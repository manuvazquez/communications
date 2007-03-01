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

#define DEBUG

PSPAlgorithm::PSPAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int firstSymbolVectorDetectedAt, double ARcoefficient): KnownChannelOrderAlgorithm(name, alphabet, L, N, K, m, channelEstimator, preamble),_inputVector(N),_stateVector(N*(m-1)),_d(smoothingLag),_startDetectionTime(preamble.cols()),_trellis(alphabet,N,m),_detectedSymbolVectors(NULL),_firstSymbolVectorDetectedAt(firstSymbolVectorDetectedAt),_ARcoefficient(ARcoefficient)
{
    if(preamble.cols() < (m-1))
        throw RuntimeException("PSPAlgorithm::PSPAlgorithm: preamble dimensions are wrong.");

    _exitStage = new PSPPath[_trellis.Nstates()];
    _arrivalStage = new PSPPath[_trellis.Nstates()];
}


PSPAlgorithm::~PSPAlgorithm()
{
    delete[] _exitStage;
    delete[] _arrivalStage;
    delete _detectedSymbolVectors;
}

void PSPAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    if(observations.cols()<(_startDetectionTime+1+_d))
        throw RuntimeException("PSPAlgorithm::Run: Not enough observations.");
}

void PSPAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
	int iProcessedObservation,iState,iBestState;

    if(observations.rows()!=_L || trainingSequence.rows()!=_N)
        throw RuntimeException("PSPAlgorithm::Run: Observations matrix or training sequence dimensions are wrong.");

    // to process the training sequence, we need both the preamble and the symbol vectors related to it
    tMatrix preambleTrainingSequence = Util::Append(_preamble,trainingSequence);


	_startDetectionTime = preambleTrainingSequence.cols();
	_detectedSymbolVectors = new tMatrix(_N,_K+_d);

    tRange rSymbolVectorsTrainingSequece(0,preambleTrainingSequence.cols()-1);

    vector<tMatrix> trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

//    	PSPPath prueba(_K+_d,0.0,preambleTrainingSequence,vector<vector<tMatrix> > (1,trainingSequenceChannelMatrices),vector<ChannelMatrixEstimator *>(1,_channelEstimator));
//    	prueba.Print();

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
    	cout << "El estado inicial es " << initialState << endl;
   	#endif

    _exitStage[initialState] = PSPPath(_K+_d,0.0,preambleTrainingSequence,vector<vector<tMatrix> > (1,trainingSequenceChannelMatrices),vector<ChannelMatrixEstimator *>(1,_channelEstimator));

	#ifdef DEBUG
    	cout << "Despues de inicializar el estado inicial" << endl;
    	cout << "El numero de estados es " << _trellis.Nstates() << endl;
   	#endif

    for(iProcessedObservation=_startDetectionTime;iProcessedObservation<_firstSymbolVectorDetectedAt;iProcessedObservation++)
    {
		#ifdef DEBUG
			cout << "-------------------- iProcessedObservation es " << iProcessedObservation << " ------------------------" << endl;
		#endif
        for(iState=0;iState<_trellis.Nstates();iState++)
        {
            if(!_exitStage[iState].IsEmpty())
            {
                DeployState(iState,observations.col(iProcessedObservation),noiseVariances[iProcessedObservation]);
            } else cout << "iState: " << iState << " estaba vacío" << endl;
        }

		#ifdef DEBUG2
			cout << "el estado 0 de _arrivalStage" << endl;
			_arrivalStage[0].Print();
		#endif

		#ifdef DEBUG
			cout << "antes de intercambiar las etapas" << endl;
		#endif

        // _arrivalStage becomes _exitStage for the next iteration
        PSPPath *aux = _exitStage;
        _exitStage = _arrivalStage;
        _arrivalStage = aux;

		#ifdef DEBUG2
			cout << "el estado 0 de _exitStage" << endl;
			_exitStage[0].Print();
		#endif

        // the _arrivalStage (old _exitStage) gets cleaned
        for(iState=0;iState<_trellis.Nstates();iState++)
			_arrivalStage[iState].Clean();

		#ifdef DEBUG2
			cout << "el estado 0 de _exitStage despues del clean" << endl;
			_exitStage[0].Print();
		#endif
    }

	#ifdef DEBUG
    	cout << "Antes de elegir el mejor estado" << initialState << endl;
   	#endif

    iBestState = BestState();

	#ifdef DEBUG
    	cout << "Despues de elegir el mejor estado" << initialState << endl;
   	#endif

    // the first detected vector is copied into "_detectedSymbolVectors"
    _detectedSymbolVectors->col(_preamble.cols()).inject(_exitStage[iBestState].GetSymbolVector(_preamble.cols()));

	#ifdef DEBUG
    	cout << "Cogiendo un vector" << initialState << endl;
   	#endif

    for( iProcessedObservation=_firstSymbolVectorDetectedAt;iProcessedObservation<_K+_d;iProcessedObservation++)
    {
        for(iState=0;iState<_trellis.Nstates();iState++)
        {
            if(!_exitStage[iState].IsEmpty())
            {
                DeployState(iState,observations.col(iProcessedObservation),noiseVariances[iProcessedObservation]);
            }
        }

        // _arrivalStage becomes _exitStage for the next iteration
        PSPPath *aux = _exitStage;
        _exitStage = _arrivalStage;
        _arrivalStage = aux;

        // the _arrivalStage (old _exitStage) gets cleaned
        for(iState=0;iState<_trellis.Nstates();iState++)
            _arrivalStage[iState].Clean();

        iBestState = BestState();

        _detectedSymbolVectors->col(iProcessedObservation-_firstSymbolVectorDetectedAt+_preamble.cols()+1).inject(_exitStage[iBestState].GetSymbolVector(iProcessedObservation-_firstSymbolVectorDetectedAt+_preamble.cols()+1));
    }

    // last detected symbol vectors are processed
    for(iProcessedObservation=_K+_d-_firstSymbolVectorDetectedAt+_preamble.cols()+1;iProcessedObservation<_K+_d;iProcessedObservation++)
        _detectedSymbolVectors->col(iProcessedObservation).inject(_exitStage[iBestState].GetSymbolVector(iProcessedObservation));
}

void PSPAlgorithm::DeployState(int iState,const tVector &observations,double noiseVariance)
{
	#ifdef DEBUG
		cout << "Al principio de DeployState (iState = " << iState << ")" << endl;
	#endif

    double newCost;
    int arrivalState;
    tVector computedObservations(_L),error(_L);

	#ifdef DEBUG
		cout << "En DeployState: antes de utilizar el estimador de canal" << endl;
	#endif

	tMatrix estimatedChannelMatrix = _exitStage[iState].GetChannelMatrixEstimator()->LastEstimatedChannelMatrix();
	estimatedChannelMatrix *= _ARcoefficient;

	#ifdef DEBUG
		cout << "En DeployState: despues de utilizar el estimador de canal" << endl;
		cout << "La matriz estimada es" << endl << estimatedChannelMatrix;
	#endif

    // "symbolVectors" will contain all the symbols involved in the current observation
    tMatrix symbolVectors(_N,_m);

	// the state determines the first "_m" symbol vectors involved in the "observations"
	_alphabet.IntToSymbolsArray(iState,_stateVector);
	for(int i=0;i<_N*(_m-1);i++)
		symbolVectors(i % _N,i / _N) = _stateVector[i];

    // now we compute the cost for each possible input
    for(int iInput=0;iInput<_trellis.NpossibleInputs();iInput++)
    {
		#ifdef DEBUG
			cout << "En DeployState: iInput es " << iInput << endl;
		#endif
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

		#ifdef DEBUG
			cout << "En DeployState: arrival state es " << arrivalState << endl;
			cout << "El nuevo coste es " << newCost << endl;
		#endif

		// if there is nothing in the arrival state
		if((_arrivalStage[arrivalState].IsEmpty()) ||
			// or there is a path whose cost is greater
			(_arrivalStage[arrivalState].GetCost() > newCost))
				// the ViterbiPath object at the arrival state is updated with that from the exit stage, the
				// new symbol vector, and the new cost
			{
				#ifdef DEBUG
					if(!_arrivalStage[arrivalState].IsEmpty())
						cout << "ECHANDO A UNO" << endl;
				#endif
				ChannelMatrixEstimator * newChannelMatrixEstimator = _exitStage[iState].GetChannelMatrixEstimator()->Clone();
				newChannelMatrixEstimator->NextMatrix(observations,symbolVectors,noiseVariance);
				_arrivalStage[arrivalState].Update(_exitStage[iState],symbolVectors.col(_m-1),newCost,vector<ChannelMatrixEstimator *>(1,newChannelMatrixEstimator));
			}
    } // for(int iInput=0;iInput<_trellis.NpossibleInputs();iInput++)


	cout << "Saliendo de DeployState" << endl;
}

void PSPAlgorithm::Process(const tMatrix &observations,vector<double> noiseVariances)
{
    // memory for the symbol vectors being detected is reserved
    _detectedSymbolVectors = new tMatrix(_N,_K+_d);
}

tMatrix PSPAlgorithm::GetDetectedSymbolVectors()
{
// 	return tMatrix(0,0);
    return (*_detectedSymbolVectors)(tRange(0,_N-1),tRange(_preamble.cols(),_K-1));
}

vector<tMatrix> PSPAlgorithm::GetEstimatedChannelMatrices()
{
	return vector<tMatrix>(0,tMatrix(0,0));
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
