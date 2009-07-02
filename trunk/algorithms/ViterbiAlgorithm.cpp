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
#include "ViterbiAlgorithm.h"

// #define DEBUG

ViterbiAlgorithm::ViterbiAlgorithm(string name, Alphabet alphabet,int L,int N, int frameLength, const StillMemoryMIMOChannel& channel,const tMatrix &preamble,int smoothingLag): KnownChannelAlgorithm(name, alphabet, L, N, frameLength,  channel),_inputVector(channel.nInputs()),_stateVector(channel.nInputs()*(channel.memory()-1)),_d(smoothingLag),_trellis(alphabet,N,channel.memory()),_preamble(preamble),_detectedSymbolVectors(NULL),rAllSymbolRows(0,_channel.nInputs()-1),rmMinus1FirstColumns(0,channel.memory()-2)
{
    if(preamble.cols() < (channel.memory()-1))
        throw RuntimeException("ViterbiAlgorithm::ViterbiAlgorithm: preamble dimensions are wrong.");

    _exitStage = new ViterbiPath[_trellis.Nstates()];
    _arrivalStage = new ViterbiPath[_trellis.Nstates()];
}


ViterbiAlgorithm::~ViterbiAlgorithm()
{
    delete[] _exitStage;
    delete[] _arrivalStage;
    delete _detectedSymbolVectors;
}

void ViterbiAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
//     Run(observations,noiseVariances,observations.cols());
    Run(observations,noiseVariances,_K+_d);
}

void ViterbiAlgorithm::Run(tMatrix observations,vector<double> noiseVariances,int firstSymbolVectorDetectedAt)
{
	const StillMemoryMIMOChannel &channel = dynamic_cast<const StillMemoryMIMOChannel &> (_channel);
    int iState,iProcessedObservation,iBestState;

    // memory for the symbol vectors being detected is reserved
    _detectedSymbolVectors = new tMatrix(channel.nInputs(),_K+_d);

    // the symbols contained in the preamble are copied into a c++ vector...
    int preambleLength = _preamble.rows()*_preamble.cols();
    vector<tSymbol> preambleVector(_nInputs*(channel.memory()-1));

    // (it must be taken into account that the number of columns of the preamble might be greater than m-1)
    int iFirstPreambleSymbolNeeded = (_preamble.cols()-(channel.memory()-1))*_nInputs;
    for(int i=iFirstPreambleSymbolNeeded;i<preambleLength;i++)
        preambleVector[i-iFirstPreambleSymbolNeeded] = _preamble(i % _preamble.rows(),i / _preamble.rows());

    // ...in order to use the method "SymbolsVectorToInt" from "Alphabet"
    int initialState = _alphabet.symbolsArray2int(preambleVector);

    _exitStage[initialState] = ViterbiPath(_K,0.0,_preamble);

    for(iProcessedObservation=_preamble.cols();iProcessedObservation<firstSymbolVectorDetectedAt;iProcessedObservation++)
    {
        for(iState=0;iState<_trellis.Nstates();iState++)
        {
            if(!_exitStage[iState].IsEmpty())
                DeployState(iState,observations.col(iProcessedObservation),channel[iProcessedObservation]);
        }

        // _arrivalStage becomes _exitStage for the next iteration
        ViterbiPath *aux = _exitStage;
        _exitStage = _arrivalStage;
        _arrivalStage = aux;

        // the _arrivalStage (old _exitStage) gets cleaned
        for(iState=0;iState<_trellis.Nstates();iState++)
			_arrivalStage[iState].Clean();
    }

    iBestState = BestState();

    // the first detected vector is copied into "_detectedSymbolVectors"
    _detectedSymbolVectors->col(_preamble.cols()).inject(_exitStage[iBestState].GetSymbolVector(_preamble.cols()));

    for( iProcessedObservation=firstSymbolVectorDetectedAt;iProcessedObservation<_K+_d;iProcessedObservation++)
    {
        for(iState=0;iState<_trellis.Nstates();iState++)
        {
            if(!_exitStage[iState].IsEmpty())
                DeployState(iState,observations.col(iProcessedObservation),channel[iProcessedObservation]);
        }

        // _arrivalStage becomes _exitStage for the next iteration
        ViterbiPath *aux = _exitStage;
        _exitStage = _arrivalStage;
        _arrivalStage = aux;

        // the _arrivalStage (old _exitStage) gets cleaned
        for(iState=0;iState<_trellis.Nstates();iState++)
            _arrivalStage[iState].Clean();

        iBestState = BestState();

        _detectedSymbolVectors->col(iProcessedObservation-firstSymbolVectorDetectedAt+_preamble.cols()+1).inject(_exitStage[iBestState].GetSymbolVector(iProcessedObservation-firstSymbolVectorDetectedAt+_preamble.cols()+1));
    }

    // last detected symbol vectors are processed
    for(iProcessedObservation=_K+_d-firstSymbolVectorDetectedAt+_preamble.cols()+1;iProcessedObservation<_K+_d;iProcessedObservation++)
        _detectedSymbolVectors->col(iProcessedObservation).inject(_exitStage[iBestState].GetSymbolVector(iProcessedObservation));
}

void ViterbiAlgorithm::DeployState(int iState,const tVector &observations,const tMatrix &channelMatrix)
{
	const StillMemoryMIMOChannel &channel = dynamic_cast<const StillMemoryMIMOChannel &> (_channel);

    double newCost;
    int arrivalState;
    tVector computedObservations(channel.nOutputs()),error(channel.nOutputs());

    // "symbolVectors" will contain all the symbols involved in the current observation
    tMatrix symbolVectors(channel.nInputs(),channel.memory());

	// the state determines the first "channel.memory()" symbol vectors involved in the "observations"
	_alphabet.int2symbolsArray(iState,_stateVector);
	for(int i=0;i<channel.nInputs()*(channel.memory()-1);i++)
		symbolVectors(i % channel.nInputs(),i / channel.nInputs()) = _stateVector[i];

    // now we compute the cost for each possible input
    for(int iInput=0;iInput<_trellis.NpossibleInputs();iInput++)
    {
        // the decimal input is converted to a symbol vector according to the alphabet
        _alphabet.int2symbolsArray(iInput,_inputVector);

        // it's copied into "symbolVectors"
        for(int i=0;i<channel.nInputs();i++)
            symbolVectors(i,channel.memory()-1) = _inputVector[i];

        // computedObservations = channelMatrix * symbolVectors(:)
        Blas_Mat_Vec_Mult(channelMatrix,Util::ToVector(symbolVectors,columnwise),computedObservations);

        // error = observations - computedObservations
        Util::add(observations,computedObservations,error,1.0,-1.0);

        newCost = _exitStage[iState].GetCost() + Blas_Dot_Prod(error,error);

        arrivalState = _trellis(iState,iInput);

		// if there is nothing in the arrival state
		if((_arrivalStage[arrivalState].IsEmpty()) ||
			// or there is a path whose cost is greater
			(_arrivalStage[arrivalState].GetCost() > newCost))
				// the ViterbiPath object at the arrival state is updated with that from the exit stage, the
				// new symbol vector, and the new cost
				_arrivalStage[arrivalState].Update(_exitStage[iState],symbolVectors.col(channel.memory()-1),newCost);
    } // for(int iInput=0;iInput<_trellis.NpossibleInputs();iInput++)

}

tMatrix ViterbiAlgorithm::getDetectedSymbolVectors()
{
    return (*_detectedSymbolVectors)(rAllSymbolRows,tRange(_preamble.cols(),_K-1));
}

void ViterbiAlgorithm::PrintStage(tStage exitOrArrival)
{
    ViterbiPath *stage;
    if(exitOrArrival == exitStage)
        stage = _exitStage;
    else
        stage = _arrivalStage;

    for(int i=0;i<_trellis.Nstates();i++)
    {
        cout << "State " << i << endl;
        if(stage[i].IsEmpty())
            cout << "Empty" << endl;
        else
		stage[i].Print();
        cout << "------------------" << endl;
    }
}
