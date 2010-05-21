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

#include "ViterbiAlgorithmWithAprioriProbabilities.h"

ViterbiAlgorithmWithAprioriProbabilities::ViterbiAlgorithmWithAprioriProbabilities(string name, Alphabet alphabet, int L, int Nr, int N, int iLastSymbolVectorToBeDetected, const StillMemoryMIMOChannel& channel, const MatrixXd& preamble, int smoothingLag, const std::vector<UsersActivityDistribution> usersActivityPdfs)
:ViterbiAlgorithm(name, alphabet,L,Nr,N, iLastSymbolVectorToBeDetected, channel,preamble,smoothingLag),
_usersActivityPdfs(usersActivityPdfs),
_extendedAlphabet(alphabet.buildNewAlphabetByAddingSymbol(0.0))
{
  if(channel.memory()>1)
	throw RuntimeException("ViterbiAlgorithmWithAprioriProbabilities::ViterbiAlgorithmWithAprioriProbabilities: algorithm is only implemented for flat channels.");
}

void ViterbiAlgorithmWithAprioriProbabilities::deployState(int iState, const VectorXd& observations, const MatrixXd& channelMatrix, const double noiseVariance)
{
    const StillMemoryMIMOChannel &channel = dynamic_cast<const StillMemoryMIMOChannel &> (_channel);

    double newCost;
    int arrivalState;

	VectorXd previousSymbolsVector = _extendedAlphabet.int2eigenVector(iState,_nInputs);
	
    // now we compute the cost for each possible input
    for(int iInput=0;iInput<_trellis->nPossibleInputs();iInput++)
    {
		// "symbolVectors" will contain the symbols involved in the current observation
		VectorXd symbolsVector = _extendedAlphabet.int2eigenVector(iInput,_nInputs);
		
        VectorXd error = observations - channelMatrix*symbolsVector;

// 		newCost = _exitStage[iState].getCost() + (error.dot(error))/(2*noiseVariance) - log(StatUtil::probXgivenY(symbolsVector,previousSymbolsVector,_usersActivityPdfs));
		newCost = _exitStage[iState].getCost() + 
				  (error.dot(error))/(2*noiseVariance) - 
				  log(StatUtil::probSymbolsVectorGivenPreviousTimeInstantUsersActivity(symbolsVector,Util::getUsersActivityFromSymbolsVector(previousSymbolsVector),_usersActivityPdfs,_alphabet.length()));

        arrivalState = (*_trellis)(iState,iInput);

        // if there is nothing in the arrival state
        if((_arrivalStage[arrivalState].isEmpty()) ||
            // or there is a path whose cost is greater
            (_arrivalStage[arrivalState].getCost() > newCost))
                // the ViterbiPath object at the arrival state is updated with that from the exit stage, the new symbol vector, and the new cost
                _arrivalStage[arrivalState].update(_exitStage[iState],symbolsVector.col(channel.memory()-1),newCost);
    } // for(int iInput=0;iInput<_trellis->nPossibleInputs();iInput++)
}

void ViterbiAlgorithmWithAprioriProbabilities::run(MatrixXd observations,vector<double> noiseVariances,int firstSymbolVectorDetectedAt)
{
    // the Trellis object is initialized (we instruct it to build a trellis assuming the memory is 2)
  _trellis = new Trellis(_extendedAlphabet,_nInputs,2);

  _exitStage = new ViterbiPath[_trellis->nStates()];
  _arrivalStage = new ViterbiPath[_trellis->nStates()];

  const StillMemoryMIMOChannel &channel = dynamic_cast<const StillMemoryMIMOChannel &> (_channel);
  
  double initialCost;
  
  // at the first time instant, a priori probabilities of the active users must be employed (instead of conditioned)
  for(int iState=0;iState<_trellis->nStates();iState++)
  {
	VectorXd symbolsVector = _extendedAlphabet.int2eigenVector(iState,_nInputs);
	
	VectorXd error = observations.col(_preamble.cols()) - channel.getTransmissionMatrix(_preamble.cols())*symbolsVector;

// 	initialCost =  (error.dot(error))/(2*noiseVariances[_preamble.cols()]) - log(StatUtil::probApriori(symbolsVector,_usersActivityPdfs));
	initialCost =  (error.dot(error))/(2*noiseVariances[_preamble.cols()]) - 
					log(StatUtil::probSymbolsVector(symbolsVector,_usersActivityPdfs,_alphabet.length()));

	_exitStage[iState] = ViterbiPath(_iLastSymbolVectorToBeDetected,initialCost,symbolsVector);
  }
  
  // first observation was already processed
  _iFirstInLoopProcessedObservation++;
  
  ViterbiAlgorithm::process(observations,noiseVariances,firstSymbolVectorDetectedAt);
}
