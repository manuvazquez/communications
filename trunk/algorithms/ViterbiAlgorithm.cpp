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

ViterbiAlgorithm::ViterbiAlgorithm(string name, Alphabet alphabet,int L,int N, int K, const StillMemoryMIMOChannel& channel,const tMatrix &preamble,int smoothingLag): KnownChannelAlgorithm(name, alphabet, L, N, K,  channel),_d(smoothingLag),_preamble(preamble),rAllSymbolRows(0,_channel.Nt()-1),rmMinus1FirstColumns(0,channel.Memory()-2)
{
    if(preamble.cols() < (channel.Memory()-1))
        throw RuntimeException("ViterbiAlgorithm::ViterbiAlgorithm: preamble dimensions are wrong.");

    _nStates = (int)pow((double)alphabet.Length(),(double)channel.Nt()*(channel.Memory()-1));
    _nPossibleInputs = (int)pow((double)alphabet.Length(),(double)channel.Nt());

    _exitStage = new tState[_nStates];
    _arrivalStage = new tState[_nStates];

    for(int i=0;i<_nStates;i++)
    {
        _exitStage[i].sequence = NULL;
        _exitStage[i].cost = 0.0;
        _arrivalStage[i].sequence = NULL;
    }

    BuildStateTransitionMatrix();
}


ViterbiAlgorithm::~ViterbiAlgorithm()
{
    for(int iState=0;iState<_nStates;iState++)
        delete[] _stateTransitionMatrix[iState];

    delete[] _stateTransitionMatrix;

    for(int i=0;i<_nStates;i++)
    {
        delete _exitStage[i].sequence;
    }

    delete[] _exitStage;
    delete[] _arrivalStage;
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
    _detectedSymbolVectors = new tMatrix(channel.Nt(),_K+_d);

    // the symbols contained in the preamble are copied into a c++ vector...
    int preambleLength = _preamble.rows()*_preamble.cols();
    vector<tSymbol> preambleVector(_N*(channel.Memory()-1));

    // (it must be taken into account that the number of columns of the preamble might be greater than m-1)
    int iFirstPreambleSymbolNeeded = (_preamble.cols()-(channel.Memory()-1))*_N;
    for(int i=iFirstPreambleSymbolNeeded;i<preambleLength;i++)
        preambleVector[i-iFirstPreambleSymbolNeeded] = _preamble(i % _preamble.rows(),i / _preamble.rows());

    // ...in order to use the method "SymbolsVectorToInt" from "Alphabet"
    int initialState = _alphabet.SymbolsArrayToInt(preambleVector);

    _exitStage[initialState].sequence = new tMatrix(_preamble);

    for( iProcessedObservation=_preamble.cols();iProcessedObservation<firstSymbolVectorDetectedAt;iProcessedObservation++)
    {
        for(iState=0;iState<_nStates;iState++)
        {
            if(_exitStage[iState].sequence!=NULL)
                DeployState(iState,observations.col(iProcessedObservation),channel[iProcessedObservation]);
        }

        // _arrivalStage becomes _exitStage for the next iteration
        tState *aux = _exitStage;
        _exitStage = _arrivalStage;
        _arrivalStage = aux;

        // the _arrivalStage (old _exitStage) gets cleaned
        for(iState=0;iState<_nStates;iState++)
        {
            delete _arrivalStage[iState].sequence;
            _arrivalStage[iState].sequence = NULL;
        }
    }

    iBestState = BestState();

    // the first detected vector is copied into "_detectedSymbolVectors"
    _detectedSymbolVectors->col(_preamble.cols()).inject((_exitStage[iBestState].sequence)->col(_preamble.cols()));

    for( iProcessedObservation=firstSymbolVectorDetectedAt;iProcessedObservation<_K+_d;iProcessedObservation++)
    {
        for(iState=0;iState<_nStates;iState++)
        {
            if(_exitStage[iState].sequence!=NULL)
                DeployState(iState,observations.col(iProcessedObservation),channel[iProcessedObservation]);
        }

        // _arrivalStage becomes _exitStage for the next iteration
        tState *aux = _exitStage;
        _exitStage = _arrivalStage;
        _arrivalStage = aux;

        // the _arrivalStage (old _exitStage) gets cleaned
        for(iState=0;iState<_nStates;iState++)
        {
            delete _arrivalStage[iState].sequence;
            _arrivalStage[iState].sequence = NULL;
        }

        iBestState = BestState();
        _detectedSymbolVectors->col(iProcessedObservation-firstSymbolVectorDetectedAt+_preamble.cols()+1).inject((_exitStage[iBestState].sequence)->col(iProcessedObservation-firstSymbolVectorDetectedAt+_preamble.cols()+1));
    }

    // last detected symbol vectors are processed
    for(iProcessedObservation=_K+_d-firstSymbolVectorDetectedAt+_preamble.cols()+1;iProcessedObservation<_K+_d;iProcessedObservation++)
        _detectedSymbolVectors->col(iProcessedObservation).inject((_exitStage[iBestState].sequence)->col(iProcessedObservation));

}

void ViterbiAlgorithm::BuildStateTransitionMatrix()
{
    _stateTransitionMatrix = new int*[_nStates];

    int alphabetLengthToTheNmMinus2 = _nStates/_nPossibleInputs;

    int iInput;
    for(int iState=0;iState<_nStates;iState++)
    {
        _stateTransitionMatrix[iState] = new int[_nPossibleInputs];
        for(iInput=0;iInput<_nPossibleInputs;iInput++)
            // computes de next state give the current one and the input (both in decimal)
            _stateTransitionMatrix[iState][iInput] = (iState % alphabetLengthToTheNmMinus2)*_nPossibleInputs + iInput;
    }
}

void ViterbiAlgorithm::DeployState(int iState,const tVector &observations,const tMatrix &channelMatrix)
{
	const StillMemoryMIMOChannel &channel = dynamic_cast<const StillMemoryMIMOChannel &> (_channel);

    double newCost;
    int arrivalState;
    tVector computedObservations(channel.Nr()),error(channel.Nr());


    int sequenceLength = (_exitStage[iState].sequence)->cols();

    // "symbolVectors" will contain all the symbols involved in the current observation
    tMatrix symbolVectors(channel.Nt(),channel.Memory());

    // the already detected symbols vectors in the sequence are copied into "symbolVectors"
    symbolVectors(rAllSymbolRows,rmMinus1FirstColumns).inject((*(_exitStage[iState].sequence))(rAllSymbolRows,tRange(sequenceLength-channel.Memory()+1,sequenceLength-1)));

    tRange rOldSymbolVectors(0,sequenceLength-1);

    // now we compute the cost for each possible input
    for(int iInput=0;iInput<_nPossibleInputs;iInput++)
    {
        // the decimal iInput is converted to a symbol vector
        vector<tSymbol> testedVector(channel.Nt());

        // according to the alphabet
        _alphabet.IntToSymbolsArray(iInput,testedVector);

        // it's copied into "symbolVectors"
        for(int i=0;i<channel.Nt();i++)
            symbolVectors(i,channel.Memory()-1) = testedVector[i];

        // computedObservations = channelMatrix * symbolVectors(:)
        Blas_Mat_Vec_Mult(channelMatrix,Util::ToVector(symbolVectors,columnwise),computedObservations);

        // error = observations - computedObservations
        Util::Add(observations,computedObservations,error,1.0,-1.0);

        newCost = _exitStage[iState].cost + Blas_Dot_Prod(error,error);

        arrivalState = _stateTransitionMatrix[iState][iInput];

        // if there is something in the arrival state
        if(_arrivalStage[arrivalState].sequence!=NULL)
        {
            // if the new cost is smaller
            if(_arrivalStage[arrivalState].cost > newCost)
            {
                // the pointer that will be overwritten is kept in order to free its memory
                tMatrix *aux = _arrivalStage[arrivalState].sequence;

                // memory for the new matrix is reserved: it will have one more column
                _arrivalStage[arrivalState].sequence = new tMatrix(channel.Nt(),sequenceLength+1);

                // the old symbol vectors are stored in the new matrix
                (*(_arrivalStage[arrivalState].sequence))(rAllSymbolRows,rOldSymbolVectors).inject(*(_exitStage[iState].sequence));

                // the new one is copied
                (_arrivalStage[arrivalState].sequence)->col(sequenceLength).inject(symbolVectors.col(channel.Memory()-1));

                // the cost is updated
                _arrivalStage[arrivalState].cost = newCost;

                // memory release
                delete aux;
            } // if(_arrivalStage[arrivalState].cost > newCost)
        } // if(_arrivalStage[arrivalState].sequence!=NULL)
        else
        {
            // memory for the new matrix is reserved: it will have one more column
            _arrivalStage[arrivalState].sequence = new tMatrix(channel.Nt(),sequenceLength+1);

            // the old symbol vectors are stored in the new matrix
            (*(_arrivalStage[arrivalState].sequence))(rAllSymbolRows,rOldSymbolVectors).inject(*(_exitStage[iState].sequence));

            // the new one is copied
            (_arrivalStage[arrivalState].sequence)->col(sequenceLength).inject(symbolVectors.col(channel.Memory()-1));

            // the cost is updated
            _arrivalStage[arrivalState].cost = newCost;
        }
    }

}

tMatrix ViterbiAlgorithm::GetDetectedSymbolVectors()
{
//     return (*_detectedSymbolVectors)(rAllSymbolRows,tRange(_channel.Memory()-1,_detectedSymbolVectors->cols()-_d-1));
    return (*_detectedSymbolVectors)(rAllSymbolRows,tRange(_preamble.cols(),_K-1));
}

void ViterbiAlgorithm::PrintStage(tStage exitOrArrival)
{
    tState *stage;
    if(exitOrArrival == exitStage)
        stage = _exitStage;
    else
        stage = _arrivalStage;

    for(int i=0;i<_nStates;i++)
    {
        cout << "State " << i << endl;
        if(stage[i].sequence==NULL)
            cout << "Empty" << endl;
        else
        {
            cout << "Sequence:" << endl << *(stage[i].sequence) << endl;
            cout << "Cost: " << stage[i].cost << endl;
        }
        cout << "------------------" << endl;
    }
}
