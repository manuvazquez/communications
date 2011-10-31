/*
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License version 2 as published by the Free Software Foundation.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public License
   along with this library; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/

#include "PSPAlgorithmWithAprioriProbabilities.h"

#include <CDMAKalmanEstimator.h>

// #define DEBUG

PSPAlgorithmWithAprioriProbabilities::PSPAlgorithmWithAprioriProbabilities(string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, uint firstSymbolVectorDetectedAt, int nSurvivors, const std::vector<UsersActivityDistribution> usersActivityPdfs):PSPAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, firstSymbolVectorDetectedAt, nSurvivors),_usersActivityPdfs(usersActivityPdfs),_extendedAlphabet(alphabet.buildNewAlphabetByAddingSymbol(0.0))
{
  if(m!=1)
	throw RuntimeException("PSPAlgorithmWithAprioriProbabilities::PSPAlgorithmWithAprioriProbabilities: this algorithm is only implemented for flat channels.");
}

void PSPAlgorithmWithAprioriProbabilities::deployState(int iState, const VectorXd& observations, double noiseVariance)
{    
    double newCost;
    int arrivalState,iDisposableSurvivor;
	
	VectorXd previousSymbolsVector = _extendedAlphabet.int2eigenVector(iState,_nInputs);

#ifdef DEBUG
	cout << "iState = " << iState << " noiseVariance = " << noiseVariance << endl;
#endif
	
    // now we compute the cost for each possible input
    for(uint iInput=0;iInput<_trellis->nPossibleInputs();iInput++)
    {
        arrivalState = (*_trellis)(iState,iInput);

		// "symbolsVector" will contain the symbols involved in the current observation
		VectorXd symbolsVector = _extendedAlphabet.int2eigenVector(iInput,_nInputs);

		for(uint iSourceSurvivor=0;iSourceSurvivor<_nSurvivors;iSourceSurvivor++)
		{
			if(_exitStage[iState][iSourceSurvivor].isEmpty())
				continue;

            VectorXd error = observations - dynamic_cast<CDMAKalmanEstimator *>(_exitStage[iState][iSourceSurvivor].getChannelMatrixEstimator())->getPredictive()*symbolsVector;
			
			newCost =  _exitStage[iState][iSourceSurvivor].getCost() + 
					  (error.dot(error))/(2*noiseVariance) - 
					  log(StatUtil::probSymbolsVectorGivenPreviousTimeInstantUsersActivity(symbolsVector,Util::getUsersActivityFromSymbolsVector(previousSymbolsVector),_usersActivityPdfs,_alphabet.length()));

            iDisposableSurvivor = disposableSurvivor(arrivalState);

			// if the given disposal survivor is empty
			if((_bestArrivingPaths[arrivalState][iDisposableSurvivor].noPathArrived()) ||
				// or its cost is greater than the computed new cost
				(_bestArrivingPaths[arrivalState][iDisposableSurvivor].getCost() > newCost))
					// the ViterbiPath object at the arrival state is updated with that from the exit stage, the
					// new symbol vector, and the new cost
				{
					_bestArrivingPaths[arrivalState][iDisposableSurvivor]._fromState = iState;
					_bestArrivingPaths[arrivalState][iDisposableSurvivor]._fromSurvivor = iSourceSurvivor;
					_bestArrivingPaths[arrivalState][iDisposableSurvivor]._input = iInput;
					_bestArrivingPaths[arrivalState][iDisposableSurvivor]._cost = newCost;
					
					// in this algorithm fields "_newSymbolVector" and "_detectedSymbolVectors" of the PSPPath are the same because
					// "_newSymbolVector" is added to the accumulated detected sequence and "_detectedSymbolVectors" is employed in
					// order to update the corresponding channel matrix estimator (Kalman estimator)
					_bestArrivingPaths[arrivalState][iDisposableSurvivor]._newSymbolVector = symbolsVector;
					_bestArrivingPaths[arrivalState][iDisposableSurvivor]._detectedSymbolVectors = symbolsVector;
				}
		} // for(int iSourceSurvivor=0;iSourceSurvivor<_nSurvivors;iSourceSurvivor++)
    } // for(int iInput=0;iInput<_trellis->nPossibleInputs();iInput++)
}

void PSPAlgorithmWithAprioriProbabilities::run(MatrixXd observations, vector< double > noiseVariances)
{
    if(observations.cols()<(_startDetectionTime+1+_d))
        throw RuntimeException("PSPAlgorithmWithAprioriProbabilities::run: not enough observations.");

    // the Trellis object is initialized
  _trellis = new Trellis(_extendedAlphabet,_nInputs,2);

  _exitStage = new PSPPath*[_trellis->nStates()];
  _arrivalStage = new PSPPath*[_trellis->nStates()];
  _bestArrivingPaths = new PSPPathCandidate*[_trellis->nStates()];

  for(uint i=0;i<_trellis->nStates();i++)
  {
	  _exitStage[i] = new PSPPath[_nSurvivors];
	  _arrivalStage[i] = new PSPPath[_nSurvivors];
	  _bestArrivingPaths[i] = new PSPPathCandidate[_nSurvivors];
  }
  
  double initialCost;
  
  // at the first time instant, a priori probabilities of the active users must be employed (instead of conditioned)
  for(uint iState=0;iState<_trellis->nStates();iState++)
  {
	VectorXd symbolsVector = _extendedAlphabet.int2eigenVector(iState,_nInputs);
	
	CDMAKalmanEstimator *clonedChannelMatrixEstimator = dynamic_cast <CDMAKalmanEstimator *> (_channelEstimator->clone());

	VectorXd error = observations.col(_startDetectionTime) - clonedChannelMatrixEstimator->lastEstimatedChannelMatrix()*symbolsVector;	

// 	initialCost =  (error.dot(error))/(2*noiseVariances[_startDetectionTime]) - log(StatUtil::probApriori(symbolsVector,_usersActivityPdfs));
	initialCost =  (error.dot(error))/(2*noiseVariances[_startDetectionTime]) - log(StatUtil::probSymbolsVector(symbolsVector,_usersActivityPdfs,_alphabet.length()));

	clonedChannelMatrixEstimator->nextMatrix(observations.col(_startDetectionTime),symbolsVector,noiseVariances[_startDetectionTime]);
	
	// only the first survivor is initialized
	_exitStage[iState][0] = PSPPath(_iLastSymbolVectorToBeDetected+_d,initialCost,symbolsVector,vector<vector<MatrixXd> > (1,vector<MatrixXd>(1,clonedChannelMatrixEstimator->lastEstimatedChannelCoefficientsMatrix())),vector<ChannelMatrixEstimator *>(1,clonedChannelMatrixEstimator));

	delete clonedChannelMatrixEstimator;
  }
  
  // first observation was already processed
  _iFirstInLoopProcessedObservation++;
	
    process(observations,noiseVariances);
}

void PSPAlgorithmWithAprioriProbabilities::run(MatrixXd observations, vector< double > noiseVariances, MatrixXd trainingSequence)
{
  throw RuntimeException("PSPAlgorithmWithAprioriProbabilities::run: this is not implemented.");
}

