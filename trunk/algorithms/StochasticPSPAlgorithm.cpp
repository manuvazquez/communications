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
#include "StochasticPSPAlgorithm.h"

// #define DEBUG

StochasticPSPAlgorithm::StochasticPSPAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int firstSymbolVectorDetectedAt, double ARcoefficient, int nSurvivors, vector<int> (*chooseSurvivors)(int,const tVector &)): KnownChannelOrderAlgorithm(name, alphabet, L, N, K, m, channelEstimator, preamble),_rAllSymbolRows(0,_N-1),_inputVector(N),_stateVector(N*(m-1)),_nSurvivors(nSurvivors),_d(smoothingLag),_startDetectionTime(preamble.cols()),_trellis(alphabet,N,m),_detectedSymbolVectors(new tMatrix(N,K+smoothingLag)),_firstSymbolVectorDetectedAt(firstSymbolVectorDetectedAt),_ARcoefficient(ARcoefficient),_chooseSurvivors(chooseSurvivors)
{
    if(preamble.cols() < (m-1))
        throw RuntimeException("StochasticPSPAlgorithm::StochasticPSPAlgorithm: preamble dimensions are wrong.");

    _exitStage = new vector<PSPPath>[_trellis.Nstates()];
    _arrivalStage = new vector<PSPPath>[_trellis.Nstates()];
    _arrivingPaths = new vector<PSPPathCandidate>[_trellis.Nstates()];

    _estimatedChannelMatrices.reserve(_K+_d-_preamble.cols());
}


StochasticPSPAlgorithm::~StochasticPSPAlgorithm()
{
    delete[] _exitStage;
    delete[] _arrivalStage;
    delete[] _arrivingPaths;
    delete _detectedSymbolVectors;
}

void StochasticPSPAlgorithm::ProcessOneObservation(const tVector &observations,double noiseVariance)
{
    #ifdef NEWDEBUG
    	cout << "-------Principio de ProcessOneObservation---------" << endl;
    #endif

	for(int iState=0;iState<_trellis.Nstates();iState++)
		// if there are survivors in this state
		if(_exitStage[iState].size()!=0)
		{
			#ifdef NEWDEBUG
				cout << "Haciendo Deploy del estado " << iState << endl;
			#endif
			DeployState(iState,observations,noiseVariance);
		}

    #ifdef NEWDEBUG
    	cout << "Hecho Deploy de todos los estados" << endl;
    #endif

	vector<double> weightsCppVector;

	for(int iState=0;iState<_trellis.Nstates();iState++)
		for(uint iSurvivor=0;iSurvivor<_arrivingPaths[iState].size();iSurvivor++)
			weightsCppVector.push_back(_arrivingPaths[iState][iSurvivor]._cost);

	// a lapackpp vector is built from the c++ one
	tVector weights(weightsCppVector.size());

	for(uint i=0;i<weightsCppVector.size();i++)
		weights(i) = weightsCppVector[i];

    #ifdef NEWDEBUG
    	cout << "Los pesos son " << endl << weights << endl;
    #endif

// 	vector<int> indexesSelectedSurvivors = StatUtil::WithoutReplacementSampling(_nSurvivors,weights);
	vector<int> indexesSelectedSurvivors = _chooseSurvivors(_nSurvivors,weights);

    #ifdef NEWDEBUG
    	cout << "Los índices seleccionados son" << endl;
    	Util::Print(indexesSelectedSurvivors);
    #endif

    #ifdef NEWDEBUG
    	cout << "Los pesos son " << endl << weights << endl;
    #endif

	double normCt = 0.0;
	for(uint i=0;i<indexesSelectedSurvivors.size();i++)
		normCt += weights(indexesSelectedSurvivors[i]);

    #ifdef NEWDEBUG
    	cout << "Calculada la cte de normalizacion: " << normCt << endl;
    #endif

	// below, the indexes are needed ordered
	sort(indexesSelectedSurvivors.begin(),indexesSelectedSurvivors.end());

    #ifdef NEWDEBUG
    	cout << "Ordenados los indices" << endl;
    #endif

	uint iSelectedSurvivor=0;
	int candidateIndex=0;
	for(int iState=0;iState<_trellis.Nstates();iState++)
	{
		vector<PSPPathCandidate>::iterator candidate = _arrivingPaths[iState].begin();
		while(candidate!=_arrivingPaths[iState].end())
		{
			#ifdef NEWDEBUG
				cout << "Estado: " << iState << " candidateIndex: " << candidateIndex << endl;
			#endif

// 			if(iSelectedSurvivor>=indexesSelectedSurvivors.size())
// 			{
// 				cout << "oooooooops" << endl;
// 				exit(1);
// 			}

			// we don't keep this survivor
			if(iSelectedSurvivor>=indexesSelectedSurvivors.size() || indexesSelectedSurvivors[iSelectedSurvivor]!=candidateIndex)
			{
				#ifdef NEWDEBUG
					cout << "Se borra el superviviente..." << endl;
				#endif
				candidate = _arrivingPaths[iState].erase(candidate);
			// we keep this survivor
			}else
			{
				#ifdef NEWDEBUG
					cout << "NO se borra el superviviente..." << endl;
				#endif
				PSPPath &sourcePath = _exitStage[candidate->_fromState][candidate->_fromSurvivor];

				#ifdef NEWDEBUG
					cout << "Localizado el Path origen (" << candidate->_fromState << "," << candidate->_fromSurvivor << ")" << endl;
				#endif

				ChannelMatrixEstimator * newChannelMatrixEstimator = sourcePath.GetChannelMatrixEstimator()->Clone();
				newChannelMatrixEstimator->NextMatrix(observations,candidate->_detectedSymbolVectors,noiseVariance);

				#ifdef NEWDEBUG
					cout << "Actualizado el estimador de canal" << endl;
					cout << "En el destino ya hay " << _arrivalStage[iState].size() << " supervivientes" << endl;
					cout << "y su capacidad es " << _arrivalStage[iState].capacity() << endl;
				#endif

				// an empty PSPPath object is created...
				_arrivalStage[iState].push_back(PSPPath());

				#ifdef NEWDEBUG
					cout << "añadido el nuevo path en el destino" << endl;
				#endif

				// ...in order to be able to call the Update method
				_arrivalStage[iState][_arrivalStage[iState].size()-1].Update(sourcePath,candidate->_newSymbolVector,candidate->_cost/normCt,vector<ChannelMatrixEstimator *>(1,newChannelMatrixEstimator));

				#ifdef NEWDEBUG
					cout << "Hecho el Update del Path recien creado" << endl;
				#endif

				iSelectedSurvivor++;
// 				// one survivor is only selected once
// 				while(indexesSelectedSurvivors[iSelectedSurvivor]==indexesSelectedSurvivors[iSelectedSurvivor-1])
// 					iSelectedSurvivor++;
				candidate++;
			}

			// anyway, the overall path candidate index must be increased
			candidateIndex++;
		}
	}

	// _arrivalStage becomes _exitStage for the next iteration
	vector<PSPPath> *aux = _exitStage;
	_exitStage = _arrivalStage;
	_arrivalStage = aux;

	// the _arrivalStage (old _exitStage) and the best arriving paths get cleaned
	for(int iState=0;iState<_trellis.Nstates();iState++)
	{
		_arrivalStage[iState].clear();
		_arrivingPaths[iState].clear();
	}

}

void StochasticPSPAlgorithm::Process(const tMatrix &observations,vector<double> noiseVariances)
{
	int iProcessedObservation,iBestState,iBestSurvivor;

    for(iProcessedObservation=_startDetectionTime;iProcessedObservation<_firstSymbolVectorDetectedAt;iProcessedObservation++)
    {
		ProcessOneObservation(observations.col(iProcessedObservation),noiseVariances[iProcessedObservation]);
    }

    BestPairStateSurvivor(iBestState,iBestSurvivor);

    // the first detected vector is copied into "_detectedSymbolVectors"...
    _detectedSymbolVectors->col(_startDetectionTime).inject(_exitStage[iBestState][iBestSurvivor].GetSymbolVector(_startDetectionTime));

    #ifdef NEWDEBUG
    	cout << "Despues del inject" << endl;
    #endif

    #ifdef NEWDEBUG
        cout << "detectado " << endl << _exitStage[iBestState][iBestSurvivor].GetSymbolVector(_startDetectionTime) << endl;
    #endif

    // ... and the first estimated channel matrix into _estimatedChannelMatrices
    _estimatedChannelMatrices.push_back(_exitStage[iBestState][iBestSurvivor].GetChannelMatrix(_startDetectionTime));
//     _estimatedChannelMatrices.push_back(_exitStage[iBestState][iBestSurvivor].GetChannelMatrixEstimator()->LastEstimatedChannelMatrix());

    #ifdef DEBUG3
    	cout << "despues de _estimatedChannelMatrices.push_back" << endl;
    #endif

    for( iProcessedObservation=_firstSymbolVectorDetectedAt;iProcessedObservation<_K+_d;iProcessedObservation++)
    {
		ProcessOneObservation(observations.col(iProcessedObservation),noiseVariances[iProcessedObservation]);

        BestPairStateSurvivor(iBestState,iBestSurvivor);

        _detectedSymbolVectors->col(iProcessedObservation-_firstSymbolVectorDetectedAt+_preamble.cols()+1).inject(_exitStage[iBestState][iBestSurvivor].GetSymbolVector(iProcessedObservation-_firstSymbolVectorDetectedAt+_preamble.cols()+1));

    	_estimatedChannelMatrices.push_back(_exitStage[iBestState][iBestSurvivor].GetChannelMatrix(iProcessedObservation));
//     	_estimatedChannelMatrices.push_back(_exitStage[iBestState][iBestSurvivor].GetChannelMatrixEstimator()->LastEstimatedChannelMatrix());
    }

    #ifdef DEBUG3
    	cout << "antes del último bucle" << endl;
    	cout << "empieza en " << _K+_d-_firstSymbolVectorDetectedAt+_startDetectionTime+1 << " y termina en " << _K+_d << endl;
    	cout << "_detectedSymbolVectors->cols(): " << _detectedSymbolVectors->cols() << endl;
    #endif

    // last detected symbol vectors are processed
    for(iProcessedObservation=_K+_d-_firstSymbolVectorDetectedAt+_startDetectionTime+1;iProcessedObservation<_K+_d;iProcessedObservation++)
    {
		#ifdef DEBUG3
			cout << "iProcessedObservation: " << iProcessedObservation << endl;
		#endif

		#ifdef DEBUG3
			cout << "eso: " << _exitStage[iBestState][iBestSurvivor]._detectedSequence->cols();
		#endif

        _detectedSymbolVectors->col(iProcessedObservation).inject(_exitStage[iBestState][iBestSurvivor].GetSymbolVector(iProcessedObservation));

		#ifdef DEBUG4
			cout << "en el medio" << iProcessedObservation << endl;
		#endif

		_estimatedChannelMatrices.push_back(_exitStage[iBestState][iBestSurvivor].GetChannelMatrix(iProcessedObservation));
// 		_estimatedChannelMatrices.push_back(_exitStage[iBestState][iBestSurvivor].GetChannelMatrixEstimator()->LastEstimatedChannelMatrix());
        #ifdef DEBUG2
            cout << "detectado " << endl << _exitStage[iBestState][iBestSurvivor].GetSymbolVector(iProcessedObservation) << endl;
        #endif
    }

    #ifdef DEBUG3
    	cout << "Final de process" << endl;
    #endif

}

void StochasticPSPAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    if(observations.cols()<(_startDetectionTime+1+_d))
        throw RuntimeException("StochasticPSPAlgorithm::Run: Not enough observations.");

    // the last N*(m-1) symbols of the preamble are copied into a c++ vector...
    int preambleLength = _preamble.rows()*_preamble.cols();
    vector<tSymbol> initialStateVector(_N*(_m-1));

    // (it must be taken into account that the number of columns of the preamble might be greater than m-1)
    int iFirstPreambleSymbolNeeded = (_preamble.cols()-(_m-1))*_N;
    for(int i=iFirstPreambleSymbolNeeded;i<preambleLength;i++)
        initialStateVector[i-iFirstPreambleSymbolNeeded] = _preamble(i % _N,i / _N);

    // ...in order to use the method "SymbolsVectorToInt" from "Alphabet" to obtain the initial state
    int initialState = _alphabet.SymbolsArrayToInt(initialStateVector);

	#ifdef NEWDEBUG
    	cout << "initialState es " << initialState << endl;
    #endif

	// the initial state is initalized
    _exitStage[initialState].push_back(PSPPath(_K+_d,1.0,_preamble,vector<vector<tMatrix> > (1,vector<tMatrix>(0)),vector<ChannelMatrixEstimator *>(1,_channelEstimator)));

	Process(observations,noiseVariances);
}

void StochasticPSPAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
	throw RuntimeException("StochasticPSPAlgorithm::Run: not implemented (with training sequence).");
}

void StochasticPSPAlgorithm::DeployState(int iState,const tVector &observations,double noiseVariance)
{
    double newCost;
    int arrivalState;
    tVector computedObservations(_L);

    // "symbolVectors" will contain all the symbols involved in the current observation
    tMatrix symbolVectors(_N,_m);

	// the state determines the first "_m" symbol vectors involved in the "observations"
	_alphabet.IntToSymbolsArray(iState,_stateVector);
	for(int i=0;i<_N*(_m-1);i++)
		symbolVectors(i % _N,i / _N) = _stateVector[i];

    // now we compute the new probability associated with each possible input
    for(int iInput=0;iInput<_trellis.NpossibleInputs();iInput++)
    {
        arrivalState = _trellis(iState,iInput);

        // the decimal input is converted to a symbol vector according to the alphabet
        _alphabet.IntToSymbolsArray(iInput,_inputVector);

        // it's copied into "symbolVectors"
        for(int i=0;i<_N;i++)
            symbolVectors(i,_m-1) = _inputVector[i];

		for(uint iSourceSurvivor=0;iSourceSurvivor<_exitStage[iState].size();iSourceSurvivor++)
		{
			tMatrix estimatedChannelMatrix = _exitStage[iState][iSourceSurvivor].GetChannelMatrixEstimator()->LastEstimatedChannelMatrix();
			estimatedChannelMatrix *= _ARcoefficient;

			// computedObservations = estimatedChannelMatrix * symbolVectors(:)
			Blas_Mat_Vec_Mult(estimatedChannelMatrix,Util::ToVector(symbolVectors,columnwise),computedObservations);

			// the cost is now the probability
			newCost = _exitStage[iState][iSourceSurvivor].GetCost()* StatUtil::NormalPdf(observations,computedObservations,noiseVariance);

			#ifdef DEBUG
				cout << "El coste (peso) actualizado es " << newCost << endl;
			#endif

			_arrivingPaths[arrivalState].push_back(PSPPathCandidate());

			// a reference to the new added PSPPathCandidate is set
			PSPPathCandidate &newSurvivor = _arrivingPaths[arrivalState][_arrivingPaths[arrivalState].size()-1];

			// the attributes are filled in
			newSurvivor._fromState = iState;
			newSurvivor._fromSurvivor = iSourceSurvivor;
			newSurvivor._input = iInput;
			newSurvivor._cost = newCost;
			newSurvivor._newSymbolVector = symbolVectors.col(_m-1);
			newSurvivor._detectedSymbolVectors = symbolVectors;

		} // for(int iSourceSurvivor=0;iSourceSurvivor<_nSurvivors;iSourceSurvivor++)
    } // for(int iInput=0;iInput<_trellis.NpossibleInputs();iInput++)

    #ifdef DEBUG
    	cout << "Una tecla..."; getchar();
    #endif
}

tMatrix StochasticPSPAlgorithm::GetDetectedSymbolVectors()
{
    return (*_detectedSymbolVectors)(_rAllSymbolRows,tRange(_preamble.cols(),_K-1));
}

vector<tMatrix> StochasticPSPAlgorithm::GetEstimatedChannelMatrices()
{
	#ifdef DEBUG3
		cout << "_d es " << _d << endl;
		cout << "_estimatedChannelMatrices: " << _estimatedChannelMatrices.size() << endl;
	#endif

	vector<tMatrix> res = _estimatedChannelMatrices;

	// the last "_d" estimated matrices need to be erased
	vector<tMatrix>::iterator tempIterator;
	tempIterator = res.end();
	tempIterator--;
	for(int i=0;i<_d;i++)
		res.erase(tempIterator);

	return res;
}

void StochasticPSPAlgorithm::BestPairStateSurvivor(int &bestState,int &bestSurvivor)
{
	if(_exitStage[0][0].IsEmpty())
		throw RuntimeException("StochasticPSPAlgorithm::BestState: [first state,first survivor] is empty.");

	double bestCost = -1.0;

	for(int iState=0;iState<_trellis.Nstates();iState++)
		for(uint iSurvivor=0;iSurvivor<_exitStage[iState].size();iSurvivor++)
		{
			#ifdef DEBUG2
				cout << "[iState: " << iState << ",iSurvivor: " << iSurvivor << "] coste = " << _exitStage[iState][iSurvivor].GetCost() << endl;
			#endif
			if(_exitStage[iState][iSurvivor].GetCost() > bestCost)
			{
				bestState = iState;
				bestSurvivor = iSurvivor;
				bestCost = _exitStage[iState][iSurvivor].GetCost();
			}
		}

	#ifdef DEBUG2
		cout << "bestState es " << bestState << " y bestSurvivor " << bestSurvivor << endl;
	#endif
}

// int StochasticPSPAlgorithm::DisposableSurvivor(int iState)
// {
//     int iWorstCost;
//     double worstCost;
//
//     if(_arrivingPaths[iState][0].NoPathArrived())
//         return 0;
//
//     iWorstCost = 0;
//     worstCost = _arrivingPaths[iState][0]._cost;
//
//     for(int iSurvivor=1;iSurvivor<_nSurvivors;iSurvivor++)
//     {
//         if(_arrivingPaths[iState][iSurvivor].NoPathArrived())
//             return iSurvivor;
//
//         if(_arrivingPaths[iState][iSurvivor]._cost > worstCost)
//         {
//             iWorstCost = iSurvivor;
//             worstCost = _arrivingPaths[iState][iSurvivor]._cost;
//         }
//     }
//
//     return iWorstCost;
// }

void StochasticPSPAlgorithm::PrintState(int iState)
{
    cout << "--------------- State " << iState << " ---------------" << endl;
    for(uint iSurvivor=0;iSurvivor<_exitStage[iState].size();iSurvivor++)
        if(!_exitStage[iState][iSurvivor].IsEmpty())
            _exitStage[iState][iSurvivor].Print();
}
