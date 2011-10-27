/***************************************************************************
 *   Copyright (C) 2006 by Manu   *
 *   manu@rustneversleeps   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "PSPAlgorithm.h"

#include <assert.h>
// #define DEBUG
// #define DEBUG_PRINT_LAST_CHANNEL_ESTIMATION

PSPAlgorithm::PSPAlgorithm(string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, uint firstSymbolVectorDetectedAt, int nSurvivors): KnownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble),_inputVector(N),_stateVector(N*(m-1)),_nSurvivors(nSurvivors),_d(smoothingLag),_startDetectionTime(preamble.cols()),_detectedSymbolVectors(new MatrixXd(N,iLastSymbolVectorToBeDetected+smoothingLag)),_firstSymbolVectorDetectedAt(firstSymbolVectorDetectedAt),_iFirstInLoopProcessedObservation(_startDetectionTime)
{
    if(preamble.cols() < (m-1))
        throw RuntimeException("PSPAlgorithm::PSPAlgorithm: preamble dimensions are wrong.");

    _estimatedChannelMatrices.reserve(_iLastSymbolVectorToBeDetected+_d-_preamble.cols());
}


PSPAlgorithm::~PSPAlgorithm()
{
    for(uint i=0;i<_trellis->nStates();i++)
    {
    	delete[] _exitStage[i];
    	delete[] _arrivalStage[i];
    	delete[] _bestArrivingPaths[i];
    }
    delete[] _exitStage;
    delete[] _arrivalStage;
    delete[] _bestArrivingPaths;
    delete _detectedSymbolVectors;
	delete _trellis;
}

void PSPAlgorithm::processOneObservation(const VectorXd &observations,double noiseVariance)
{
	uint iState,iSurvivor;

	for(iState=0;iState<_trellis->nStates();iState++)
	{
		// if the first survivor is empty, we assume that so are the remaining ones
		if(!_exitStage[iState][0].isEmpty())
			deployState(iState,observations,noiseVariance);
	}

	// the best paths arriving at each state are generated from the stored PathCandidate
	for(iState=0;iState<_trellis->nStates();iState++)
	{
		for(iSurvivor=0;iSurvivor<_nSurvivors;iSurvivor++)
		{
			PSPPathCandidate &bestPathCandidate = _bestArrivingPaths[iState][iSurvivor];
			
			if(bestPathCandidate.noPathArrived())
				continue;

			PSPPath &sourcePath = _exitStage[bestPathCandidate._fromState][bestPathCandidate._fromSurvivor];
			ChannelMatrixEstimator * newChannelMatrixEstimator = sourcePath.getChannelMatrixEstimator()->clone();

            newChannelMatrixEstimator->nextMatrix(observations,bestPathCandidate._detectedSymbolVectors,noiseVariance);
			
#ifdef DEBUG_PRINT_UPDATED_CHANNEL_MATRICES
			cout << "updated channel matrix estimator" << endl << newChannelMatrixEstimator->lastEstimatedChannelMatrix_eigen() << endl;
#endif

			_arrivalStage[iState][iSurvivor].update(sourcePath,bestPathCandidate._newSymbolVector,bestPathCandidate._cost,vector<ChannelMatrixEstimator *>(1,newChannelMatrixEstimator));
		}
	}

	// _arrivalStage becomes _exitStage for the next iteration
	PSPPath **aux = _exitStage;
	_exitStage = _arrivalStage;
	_arrivalStage = aux;

	// the _arrivalStage (old _exitStage) and the best arriving paths get cleaned
	for(iState=0;iState<_trellis->nStates();iState++)
	{
		for(iSurvivor=0;iSurvivor<_nSurvivors;iSurvivor++)
		{
			_arrivalStage[iState][iSurvivor].clean();
			_bestArrivingPaths[iState][iSurvivor].clean();
		}
	}
}

void PSPAlgorithm::process(const MatrixXd &observations,vector<double> noiseVariances)
{
	uint iProcessedObservation,iBestState,iBestSurvivor;

//     for(iProcessedObservation=_startDetectionTime;iProcessedObservation<_firstSymbolVectorDetectedAt;iProcessedObservation++)
    for(iProcessedObservation=_iFirstInLoopProcessedObservation;iProcessedObservation<_firstSymbolVectorDetectedAt;iProcessedObservation++)
	{
		processOneObservation(observations.col(iProcessedObservation),noiseVariances[iProcessedObservation]);
#ifdef DEBUG
		cout << "============== iProcessedObservation = " << iProcessedObservation << " =====================" << endl;
		getchar();
#endif
	}

    bestPairStateSurvivor(iBestState,iBestSurvivor);

    // the first detected vector is copied into "_detectedSymbolVectors"...
    _detectedSymbolVectors->col(_startDetectionTime) = _exitStage[iBestState][iBestSurvivor].getSymbolVector(_startDetectionTime);

    // ... and the first estimated channel matrix into _estimatedChannelMatrices
    _estimatedChannelMatrices.push_back(_exitStage[iBestState][iBestSurvivor].getChannelMatrix(_startDetectionTime));

    for(iProcessedObservation=_firstSymbolVectorDetectedAt;iProcessedObservation<_iLastSymbolVectorToBeDetected+_d;iProcessedObservation++)
    {
		processOneObservation(observations.col(iProcessedObservation),noiseVariances[iProcessedObservation]);

        bestPairStateSurvivor(iBestState,iBestSurvivor);

        _detectedSymbolVectors->col(iProcessedObservation-_firstSymbolVectorDetectedAt+_preamble.cols()+1) = _exitStage[iBestState][iBestSurvivor].getSymbolVector(iProcessedObservation-_firstSymbolVectorDetectedAt+_preamble.cols()+1);


    	_estimatedChannelMatrices.push_back(_exitStage[iBestState][iBestSurvivor].getChannelMatrix(iProcessedObservation));
    }

#ifdef DEBUG_PRINT_LAST_CHANNEL_ESTIMATION
	cout << "last estimated channel matrix is" << endl << _exitStage[iBestState][iBestSurvivor].getChannelMatrixEstimator()->lastEstimatedChannelMatrix() << endl;
#endif

    // last detected symbol vectors are processed
    for(iProcessedObservation=_iLastSymbolVectorToBeDetected+_d-_firstSymbolVectorDetectedAt+_startDetectionTime+1;iProcessedObservation<_iLastSymbolVectorToBeDetected+_d;iProcessedObservation++)
    {
        _detectedSymbolVectors->col(iProcessedObservation) = _exitStage[iBestState][iBestSurvivor].getSymbolVector(iProcessedObservation);

		_estimatedChannelMatrices.push_back(_exitStage[iBestState][iBestSurvivor].getChannelMatrix(iProcessedObservation));
    }
}

void PSPAlgorithm::run(MatrixXd observations,vector<double> noiseVariances)
{
    if(observations.cols()<(_startDetectionTime+1+_d))
        throw RuntimeException("PSPAlgorithm::Run: not enough observations.");

	initializeTrellis();

    // the last N*(m-1) symbols of the preamble are copied into a c++ vector...
    uint preambleLength = _preamble.rows()*_preamble.cols();
    vector<tSymbol> initialStateVector(_nInputs*(_channelOrder-1));

    // (it must be taken into account that the number of columns of the preamble might be greater than m-1)
	assert(_preamble.cols()-(_channelOrder-1)>=0);
    uint iFirstPreambleSymbolNeeded = (_preamble.cols()-(_channelOrder-1))*_nInputs;
    for(uint i=iFirstPreambleSymbolNeeded;i<preambleLength;i++)
        initialStateVector[i-iFirstPreambleSymbolNeeded] = _preamble(i % _nInputs,i / _nInputs);

    // ...in order to use the method "SymbolsVectorToInt" from "Alphabet" to obtain the initial state
    int initialState = _alphabet.symbolsArray2int(initialStateVector);

    // the initial state is initalized
    _exitStage[initialState][0] = PSPPath(_iLastSymbolVectorToBeDetected+_d,0.0,_preamble,vector<vector<MatrixXd> > (1,vector<MatrixXd>(0)),vector<ChannelMatrixEstimator *>(1,_channelEstimator));

    process(observations,noiseVariances);
}

void PSPAlgorithm::run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence)
{
    if(observations.rows()!=_nOutputs || trainingSequence.rows()!=_nInputs)
        throw RuntimeException("PSPAlgorithm::Run: Observations matrix or training sequence dimensions are wrong.");

	initializeTrellis();

    // to process the training sequence, we need both the preamble and the symbol vectors related to it
    MatrixXd preambleTrainingSequence(trainingSequence.rows(),_preamble.cols()+trainingSequence.cols());
    preambleTrainingSequence << _preamble,trainingSequence;

    _startDetectionTime = preambleTrainingSequence.cols();
	
	// this is needed for PSPAlgorithmWithAprioriProbabilities
	_iFirstInLoopProcessedObservation = _startDetectionTime;

    vector<MatrixXd> trainingSequenceChannelMatrices = _channelEstimator->nextMatricesFromObservationsSequence(observations,noiseVariances,preambleTrainingSequence,_preamble.cols(),_startDetectionTime);

    // known symbol vectors are copied into the the vector with the final detected ones
    _detectedSymbolVectors->block(0,_preamble.cols(),_nInputs,_startDetectionTime-_preamble.cols()) = trainingSequence;

    // and so the channel matrices
    for(uint i=0;i<trainingSequenceChannelMatrices.size();i++)
        _estimatedChannelMatrices.push_back(trainingSequenceChannelMatrices[i]);

    // the last N*(m-1) symbols of the training sequence are copied into a c++ vector...
    int preambleTrainingSequenceLength = preambleTrainingSequence.size();
    vector<tSymbol> initialStateVector(_nInputs*(_channelOrder-1));

    // (it must be taken into account that the number of columns of the preamble might be greater than m-1)
    int iFirstPreambleSymbolNeeded = (preambleTrainingSequence.cols()-(_channelOrder-1))*_nInputs;
    for(int i=iFirstPreambleSymbolNeeded;i<preambleTrainingSequenceLength;i++)
        initialStateVector[i-iFirstPreambleSymbolNeeded] = preambleTrainingSequence(i % _nInputs,i / _nInputs);

    // ...in order to use the method "SymbolsVectorToInt" from "Alphabet" to obtain the initial state
    int initialState = _alphabet.symbolsArray2int(initialStateVector);

    // the initial state is initalized
    _exitStage[initialState][0] = PSPPath(_iLastSymbolVectorToBeDetected+_d,0.0,preambleTrainingSequence,vector<vector<MatrixXd> > (1,trainingSequenceChannelMatrices),vector<ChannelMatrixEstimator *>(1,_channelEstimator));

    process(observations,noiseVariances);
}

void PSPAlgorithm::deployState(int iState,const VectorXd &observations,double noiseVariance)
{
    double newCost;
    int arrivalState,iDisposableSurvivor;

    // "symbolVectors" will contain all the symbols involved in the current observation
    MatrixXd symbolVectors(_nInputs,_channelOrder);

	// the state determines the first "_channelOrder"-1 symbol vectors involved in the "observations"
	_alphabet.int2symbolsArray(iState,_stateVector);
	for(uint i=0;i<_nInputs*(_channelOrder-1);i++)
		symbolVectors(i % _nInputs,i / _nInputs) = _stateVector[i];

    // now we compute the cost for each possible input
    for(uint iInput=0;iInput<_trellis->nPossibleInputs();iInput++)
    {
        arrivalState = (*_trellis)(iState,iInput);

        // the decimal input is converted to a symbol vector according to the alphabet
        _alphabet.int2symbolsArray(iInput,_inputVector);

        // it's copied into "symbolVectors"
        for(uint i=0;i<_nInputs;i++)
            symbolVectors(i,_channelOrder-1) = _inputVector[i];

		for(uint iSourceSurvivor=0;iSourceSurvivor<_nSurvivors;iSourceSurvivor++)
		{
			if(_exitStage[iState][iSourceSurvivor].isEmpty())
				continue;

//             VectorXd error = observations - _ARcoefficient*_exitStage[iState][iSourceSurvivor].getChannelMatrixEstimator()->lastEstimatedChannelMatrix_eigen()*Util::toVector(symbolVectors,columnwise);

            VectorXd error = observations - dynamic_cast<KalmanEstimator *>(_exitStage[iState][iSourceSurvivor].getChannelMatrixEstimator())->getPredictive()*Util::toVector(symbolVectors,columnwise);

#ifdef DEBUG
			cout << "las estimated chanel matrix" << endl << dynamic_cast<KalmanEstimator *>(_exitStage[iState][iSourceSurvivor].getChannelMatrixEstimator())->getPredictive() << endl;
#endif

            newCost = _exitStage[iState][iSourceSurvivor].getCost() + error.dot(error);

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
					_bestArrivingPaths[arrivalState][iDisposableSurvivor]._newSymbolVector = symbolVectors.col(_channelOrder-1);
					_bestArrivingPaths[arrivalState][iDisposableSurvivor]._detectedSymbolVectors = symbolVectors;
				}
		} // for(int iSourceSurvivor=0;iSourceSurvivor<_nSurvivors;iSourceSurvivor++)
    } // for(int iInput=0;iInput<_trellis->nPossibleInputs();iInput++)
}

MatrixXd PSPAlgorithm::getDetectedSymbolVectors()
{
    return _detectedSymbolVectors->block(0,_preamble.cols(),_nInputs,_iLastSymbolVectorToBeDetected-_preamble.cols());
}

vector<MatrixXd> PSPAlgorithm::getEstimatedChannelMatrices()
{
    vector<MatrixXd> res = _estimatedChannelMatrices;

    // the last "_d" estimated matrices need to be erased
    vector<MatrixXd>::iterator tempIterator;
    tempIterator = res.end();
    tempIterator--;
    for(uint i=0;i<_d;i++)
        res.erase(tempIterator);

    return res;
}

void PSPAlgorithm::bestPairStateSurvivor(uint& bestState, uint& bestSurvivor)
{
	if(_exitStage[0][0].isEmpty())
		throw RuntimeException("PSPAlgorithm::bestState: [first state,first survivor] is empty.");

	uint iSurvivor;
	bestState = 0;
	bestSurvivor = 0;
	double bestCost = _exitStage[bestState][bestSurvivor].getCost();

	for(uint iState=1;iState<_trellis->nStates();iState++)
		for(iSurvivor=0;iSurvivor<_nSurvivors;iSurvivor++)
		{
			if(_exitStage[iState][iSurvivor].getCost() < bestCost)
			{
				bestState = iState;
				bestSurvivor = iSurvivor;
				bestCost = _exitStage[iState][iSurvivor].getCost();
			}
		}
}

int PSPAlgorithm::disposableSurvivor(int iState)
{
    int iWorstCost;
    double worstCost;

	// if there are no survivors yet
    if(_bestArrivingPaths[iState][0].noPathArrived())
		// the first slot is ok
        return 0;

    iWorstCost = 0;
    worstCost = _bestArrivingPaths[iState][0]._cost;

    for(uint iSurvivor=1;iSurvivor<_nSurvivors;iSurvivor++)
    {
        if(_bestArrivingPaths[iState][iSurvivor].noPathArrived())
            return iSurvivor;

        if(_bestArrivingPaths[iState][iSurvivor]._cost > worstCost)
        {
            iWorstCost = iSurvivor;
            worstCost = _bestArrivingPaths[iState][iSurvivor]._cost;
        }
    }

    return iWorstCost;
}

void PSPAlgorithm::printState(int iState)
{
    cout << "--------------- State " << iState << " ---------------" << endl;
    for(uint iSurvivor=0;iSurvivor<_nSurvivors;iSurvivor++)
        if(!_exitStage[iState][iSurvivor].isEmpty())
            _exitStage[iState][iSurvivor].print();
}

void PSPAlgorithm::initializeTrellis()
{
  _trellis = new Trellis(_alphabet,_nInputs,_channelOrder);
  
  _exitStage = new PSPPath*[_trellis->nStates()];
  _arrivalStage = new PSPPath*[_trellis->nStates()];
  _bestArrivingPaths = new PSPPathCandidate*[_trellis->nStates()];

  for(uint i=0;i<_trellis->nStates();i++)
  {
	  _exitStage[i] = new PSPPath[_nSurvivors];
	  _arrivalStage[i] = new PSPPath[_nSurvivors];
	  _bestArrivingPaths[i] = new PSPPathCandidate[_nSurvivors];
  }
}