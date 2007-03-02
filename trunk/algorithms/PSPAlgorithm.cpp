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

// #define DEBUG

PSPAlgorithm::PSPAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int firstSymbolVectorDetectedAt, double ARcoefficient): KnownChannelOrderAlgorithm(name, alphabet, L, N, K, m, channelEstimator, preamble),_rAllSymbolRows(0,_N-1),_inputVector(N),_stateVector(N*(m-1)),_d(smoothingLag),_startDetectionTime(preamble.cols()),_trellis(alphabet,N,m),_detectedSymbolVectors(new tMatrix(N,K+smoothingLag)),_firstSymbolVectorDetectedAt(firstSymbolVectorDetectedAt),_ARcoefficient(ARcoefficient)
{
    if(preamble.cols() < (m-1))
        throw RuntimeException("PSPAlgorithm::PSPAlgorithm: preamble dimensions are wrong.");

    _exitStage = new PSPPath[_trellis.Nstates()];
    _arrivalStage = new PSPPath[_trellis.Nstates()];
    _bestArrivingPaths = new PSPPathCandidate[_trellis.Nstates()];

    _estimatedChannelMatrices.reserve(_K+_d-_preamble.cols());
}


PSPAlgorithm::~PSPAlgorithm()
{
    delete[] _exitStage;
    delete[] _arrivalStage;
    delete[] _bestArrivingPaths;
    delete _detectedSymbolVectors;
}

void PSPAlgorithm::ProcessOneObservation(const tVector &observations,double noiseVariance)
{
	int iState;

	for(iState=0;iState<_trellis.Nstates();iState++)
	{
		if(!_exitStage[iState].IsEmpty())
		{
			DeployState(iState,observations,noiseVariance);
		}
	}

	// the best path arriving at each state is generated from the stored PathSourceId
	for(iState=0;iState<_trellis.Nstates();iState++)
	{
		if(_bestArrivingPaths[iState].NoPathArrived())
			continue;

		ChannelMatrixEstimator * newChannelMatrixEstimator = _exitStage[_bestArrivingPaths[iState]._fromState].GetChannelMatrixEstimator()->Clone();

		newChannelMatrixEstimator->NextMatrix(observations,_bestArrivingPaths[iState]._detectedSymbolVectors,noiseVariance);
		_arrivalStage[iState].Update(_exitStage[_bestArrivingPaths[iState]._fromState],_bestArrivingPaths[iState]._newSymbolVector,_bestArrivingPaths[iState]._cost,vector<ChannelMatrixEstimator *>(1,newChannelMatrixEstimator));
	}

	// _arrivalStage becomes _exitStage for the next iteration
	PSPPath *aux = _exitStage;
	_exitStage = _arrivalStage;
	_arrivalStage = aux;

	// the _arrivalStage (old _exitStage) and the best arriving paths get cleaned
	for(iState=0;iState<_trellis.Nstates();iState++)
	{
		_arrivalStage[iState].Clean();
		_bestArrivingPaths[iState].Clean();
	}
}

void PSPAlgorithm::Process(const tMatrix &observations,vector<double> noiseVariances)
{
	int iProcessedObservation,iBestState;

    for(iProcessedObservation=_startDetectionTime;iProcessedObservation<_firstSymbolVectorDetectedAt;iProcessedObservation++)
    {
//     	cout << "iProcessedObservation = " << iProcessedObservation << endl;
		ProcessOneObservation(observations.col(iProcessedObservation),noiseVariances[iProcessedObservation]);
    } // for(iProcessedObservation=_startDetectionTime;iProcessedObservation<_firstSymbolVectorDetectedAt;iProcessedObservation++)

    iBestState = BestState();

    // the first detected vector is copied into "_detectedSymbolVectors"...
    _detectedSymbolVectors->col(_startDetectionTime).inject(_exitStage[iBestState].GetSymbolVector(_startDetectionTime));

    // ... and the first estimated channel matrix into _estimatedChannelMatrices
    _estimatedChannelMatrices.push_back(_exitStage[iBestState].GetChannelMatrixEstimator()->LastEstimatedChannelMatrix());

    for( iProcessedObservation=_firstSymbolVectorDetectedAt;iProcessedObservation<_K+_d;iProcessedObservation++)
    {
		ProcessOneObservation(observations.col(iProcessedObservation),noiseVariances[iProcessedObservation]);

        iBestState = BestState();

        _detectedSymbolVectors->col(iProcessedObservation-_firstSymbolVectorDetectedAt+_preamble.cols()+1).inject(_exitStage[iBestState].GetSymbolVector(iProcessedObservation-_firstSymbolVectorDetectedAt+_preamble.cols()+1));
    	_estimatedChannelMatrices.push_back(_exitStage[iBestState].GetChannelMatrixEstimator()->LastEstimatedChannelMatrix());
    }

    // last detected symbol vectors are processed
    for(iProcessedObservation=_K+_d-_firstSymbolVectorDetectedAt+_startDetectionTime+1;iProcessedObservation<_K+_d;iProcessedObservation++)
    {
        _detectedSymbolVectors->col(iProcessedObservation).inject(_exitStage[iBestState].GetSymbolVector(iProcessedObservation));
		_estimatedChannelMatrices.push_back(_exitStage[iBestState].GetChannelMatrixEstimator()->LastEstimatedChannelMatrix());
    }
}

void PSPAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    if(observations.cols()<(_startDetectionTime+1+_d))
        throw RuntimeException("PSPAlgorithm::Run: Not enough observations.");
}

void PSPAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
    if(observations.rows()!=_L || trainingSequence.rows()!=_N)
        throw RuntimeException("PSPAlgorithm::Run: Observations matrix or training sequence dimensions are wrong.");

    // to process the training sequence, we need both the preamble and the symbol vectors related to it
    tMatrix preambleTrainingSequence = Util::Append(_preamble,trainingSequence);

	_startDetectionTime = preambleTrainingSequence.cols();

    vector<tMatrix> trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

	// known symbol vectors are copied into the the vector with the final detected ones
	(*_detectedSymbolVectors)(_rAllSymbolRows,tRange(_preamble.cols(),_startDetectionTime-1)).inject(trainingSequence);

	// and so the channel matrices
	for(uint i=_preamble.cols();i<trainingSequenceChannelMatrices.size();i++)
		_estimatedChannelMatrices.push_back(trainingSequenceChannelMatrices[i]);

    // the symbols contained in the preamble are copied into a c++ vector...
    int preambleTrainingSequenceLength = preambleTrainingSequence.rows()*preambleTrainingSequence.cols();
    vector<tSymbol> initialStateVector(_N*(_m-1));

    // (it must be taken into account that the number of columns of the preamble might be greater than m-1)
    int iFirstPreambleSymbolNeeded = (preambleTrainingSequence.cols()-(_m-1))*_N;
    for(int i=iFirstPreambleSymbolNeeded;i<preambleTrainingSequenceLength;i++)
        initialStateVector[i-iFirstPreambleSymbolNeeded] = preambleTrainingSequence(i % _N,i / _N);

    // ...in order to use the method "SymbolsVectorToInt" from "Alphabet"
    int initialState = _alphabet.SymbolsArrayToInt(initialStateVector);

	#ifdef DEBUG
    	cout << "initialState es " << initialState << endl;
    #endif

	// the initial state is initalized
    _exitStage[initialState] = PSPPath(_K+_d,0.0,preambleTrainingSequence,vector<vector<tMatrix> > (1,trainingSequenceChannelMatrices),vector<ChannelMatrixEstimator *>(1,_channelEstimator));

	Process(observations,noiseVariances);
}

void PSPAlgorithm::DeployState(int iState,const tVector &observations,double noiseVariance)
{
    double newCost;
    int arrivalState;
    tVector computedObservations(_L),error(_L);

	tMatrix estimatedChannelMatrix = _exitStage[iState].GetChannelMatrixEstimator()->LastEstimatedChannelMatrix();
	estimatedChannelMatrix *= _ARcoefficient;

    // "symbolVectors" will contain all the symbols involved in the current observation
    tMatrix symbolVectors(_N,_m);

	// the state determines the first "_m" symbol vectors involved in the "observations"
	_alphabet.IntToSymbolsArray(iState,_stateVector);
	for(int i=0;i<_N*(_m-1);i++)
		symbolVectors(i % _N,i / _N) = _stateVector[i];

    // now we compute the cost for each possible input
    for(int iInput=0;iInput<_trellis.NpossibleInputs();iInput++)
    {
        // the decimal input is converted to a symbol vector according to the alphabet
        _alphabet.IntToSymbolsArray(iInput,_inputVector);

        // it's copied into "symbolVectors"
        for(int i=0;i<_N;i++)
            symbolVectors(i,_m-1) = _inputVector[i];

        // computedObservations = estimatedChannelMatrix * symbolVectors(:)
        Blas_Mat_Vec_Mult(estimatedChannelMatrix,Util::ToVector(symbolVectors,columnwise),computedObservations);

        // error = observations - computedObservations
        Util::Add(observations,computedObservations,error,1.0,-1.0);

        newCost = _exitStage[iState].GetCost() + Blas_Dot_Prod(error,error);

        arrivalState = _trellis(iState,iInput);

		// if there is nothing in the arrival state
		if((_bestArrivingPaths[arrivalState].NoPathArrived()) ||
			// or there is a path whose cost is greater
			(_bestArrivingPaths[arrivalState].GetCost() > newCost))
				// the ViterbiPath object at the arrival state is updated with that from the exit stage, the
				// new symbol vector, and the new cost
			{
				_bestArrivingPaths[arrivalState]._fromState = iState;
				_bestArrivingPaths[arrivalState]._input = iInput;
				_bestArrivingPaths[arrivalState]._cost = newCost;
				_bestArrivingPaths[arrivalState]._newSymbolVector = symbolVectors.col(_m-1);
				_bestArrivingPaths[arrivalState]._detectedSymbolVectors = symbolVectors;
			}
    } // for(int iInput=0;iInput<_trellis.NpossibleInputs();iInput++)
}

tMatrix PSPAlgorithm::GetDetectedSymbolVectors()
{
    return (*_detectedSymbolVectors)(_rAllSymbolRows,tRange(_preamble.cols(),_K-1));
}

vector<tMatrix> PSPAlgorithm::GetEstimatedChannelMatrices()
{
	// the last "_d" estimated matrices need to be erased
	vector<tMatrix>::iterator tempIterator;
	tempIterator = _estimatedChannelMatrices.end();
	tempIterator--;
	for(int i=0;i<_d;i++)
		_estimatedChannelMatrices.erase(tempIterator);

	return _estimatedChannelMatrices;
}

int PSPAlgorithm::BestState()
{
	if(_exitStage[0].IsEmpty())
		throw RuntimeException("PSPAlgorithm::BestState: first state is empty.");

	int bestState = 0;
	double bestCost = _exitStage[0].GetCost();

	for(int iState=1;iState<_trellis.Nstates();iState++)
		if(_exitStage[iState].GetCost() < bestCost)
		{
			bestState = iState;
			bestCost = _exitStage[iState].GetCost();
		}
	return bestState;
}
