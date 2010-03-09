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

ViterbiAlgorithm::ViterbiAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected, const StillMemoryMIMOChannel& channel,const MatrixXd &preamble,int smoothingLag): KnownChannelAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected,  channel),_d(smoothingLag),_trellis(NULL),_preamble(preamble),_detectedSymbolVectors(NULL),_iFirstInLoopProcessedObservation(_preamble.cols())
{
    if(preamble.cols() < (channel.memory()-1))
        throw RuntimeException("ViterbiAlgorithm::ViterbiAlgorithm: preamble dimensions are wrong.");
}


ViterbiAlgorithm::~ViterbiAlgorithm()
{
    delete[] _exitStage;
    delete[] _arrivalStage;
    delete _detectedSymbolVectors;
	delete _trellis;
}

void ViterbiAlgorithm::run(MatrixXd observations,vector<double> noiseVariances)
{
    run(observations,noiseVariances,_iLastSymbolVectorToBeDetected+_d);
}

void ViterbiAlgorithm::run(MatrixXd observations,vector<double> noiseVariances,int firstSymbolVectorDetectedAt)
{
  const StillMemoryMIMOChannel &channel = dynamic_cast<const StillMemoryMIMOChannel &> (_channel);
  
  if(channel.memory()==1)
	throw RuntimeException("ViterbiAlgorithm::initializeTrellis: algorithm is only implemented for frequency-selective channels.");
  
  // the Trellis object is initialized
  _trellis = new Trellis(_alphabet,_nInputs,channel.memory());

  _exitStage = new ViterbiPath[_trellis->Nstates()];
  _arrivalStage = new ViterbiPath[_trellis->Nstates()];

    // the symbols contained in the preamble are copied into a c++ vector...
    vector<tSymbol> preambleVector(_nInputs*(channel.memory()-1));

    // (it must be taken into account that the number of columns of the preamble might be greater than m-1)
    int iFirstPreambleSymbolNeeded = (_preamble.cols()-(channel.memory()-1))*_nInputs;
    for(int i=iFirstPreambleSymbolNeeded;i<_preamble.size();i++)
        preambleVector[i-iFirstPreambleSymbolNeeded] = _preamble(i % _preamble.rows(),i / _preamble.rows());
  
	// the sequence of symbols is converted to a number using the corresponding method from "Alphabet"
    int initialState = _alphabet.symbolsArray2int(preambleVector);

    _exitStage[initialState] = ViterbiPath(_iLastSymbolVectorToBeDetected,0.0,_preamble);
  
  process(observations,noiseVariances,firstSymbolVectorDetectedAt);
}

void ViterbiAlgorithm::process(MatrixXd observations,vector<double> noiseVariances,int firstSymbolVectorDetectedAt)
{  
    const StillMemoryMIMOChannel &channel = dynamic_cast<const StillMemoryMIMOChannel &> (_channel);
    int iState,iProcessedObservation,iBestState;

    // memory for the symbol vectors being detected is reserved
    _detectedSymbolVectors = new MatrixXd(channel.nInputs(),_iLastSymbolVectorToBeDetected+_d);

    for(iProcessedObservation=_iFirstInLoopProcessedObservation;iProcessedObservation<firstSymbolVectorDetectedAt;iProcessedObservation++)
    {
        for(iState=0;iState<_trellis->Nstates();iState++)
        {
            if(!_exitStage[iState].IsEmpty())
                DeployState(iState,observations.col(iProcessedObservation),channel.getTransmissionMatrix(iProcessedObservation),noiseVariances[iProcessedObservation]);
        }

        // _arrivalStage becomes _exitStage for the next iteration
		swapStages();
    }

    iBestState = BestState();

    // the first detected vector is copied into "_detectedSymbolVectors"
    _detectedSymbolVectors->col(_preamble.cols()) = _exitStage[iBestState].getSymbolVector(_preamble.cols());

    for( iProcessedObservation=firstSymbolVectorDetectedAt;iProcessedObservation<_iLastSymbolVectorToBeDetected+_d;iProcessedObservation++)
    {
        for(iState=0;iState<_trellis->Nstates();iState++)
        {
            if(!_exitStage[iState].IsEmpty())
                DeployState(iState,observations.col(iProcessedObservation),channel.getTransmissionMatrix(iProcessedObservation),noiseVariances[iProcessedObservation]);
        }

        // _arrivalStage becomes _exitStage for the next iteration
		swapStages();

        iBestState = BestState();

        _detectedSymbolVectors->col(iProcessedObservation-firstSymbolVectorDetectedAt+_preamble.cols()+1) = _exitStage[iBestState].getSymbolVector(iProcessedObservation-firstSymbolVectorDetectedAt+_preamble.cols()+1);
    }

#ifdef DEBUG
		cout << "_detectedSymbolVectors = " << endl << *_detectedSymbolVectors << endl;
		cout << "_exitStage[iBestState]" << endl << _exitStage[iBestState].getDetectedSequence() << endl;
#endif

    // last detected symbol vectors are processed
    for(iProcessedObservation=_iLastSymbolVectorToBeDetected+_d-firstSymbolVectorDetectedAt+_preamble.cols()+1;iProcessedObservation<_iLastSymbolVectorToBeDetected+_d;iProcessedObservation++)
	{
#ifdef DEBUG
		cout << "iProcessedObservation = " << iProcessedObservation << endl;
		cout << "vector to add" << endl << _exitStage[iBestState].getSymbolVector(iProcessedObservation) << endl;
#endif
        _detectedSymbolVectors->col(iProcessedObservation) = _exitStage[iBestState].getSymbolVector(iProcessedObservation);
	}
}

void ViterbiAlgorithm::DeployState(int iState,const VectorXd &observations,const MatrixXd &channelMatrix,const double noiseVariance)
{
    const StillMemoryMIMOChannel &channel = dynamic_cast<const StillMemoryMIMOChannel &> (_channel);

    double newCost;
    int arrivalState;

    // "symbolVectors" will contain all the symbols involved in the current observation
    MatrixXd symbolVectors(channel.nInputs(),channel.memory());

	// the state determines the first "channel.memory()" - 1 symbol vectors involved in the "observations"
	symbolVectors.block(0,0,channel.nInputs(),channel.memory()-1) = _alphabet.int2eigenMatrix(iState,channel.nInputs(),channel.memory()-1);
	
    // now we compute the cost for each possible input
    for(int iInput=0;iInput<_trellis->NpossibleInputs();iInput++)
    {
		symbolVectors.col(channel.memory()-1) = _alphabet.int2eigenVector(iInput,channel.nInputs());
		
        VectorXd error = observations - channelMatrix*Util::toVector(symbolVectors,columnwise);

        newCost = _exitStage[iState].getCost() + error.dot(error);

        arrivalState = (*_trellis)(iState,iInput);

        // if there is nothing in the arrival state
        if((_arrivalStage[arrivalState].IsEmpty()) ||
            // or there is a path whose cost is greater
            (_arrivalStage[arrivalState].getCost() > newCost))
                // the ViterbiPath object at the arrival state is updated with that from the exit stage, the new symbol vector, and the new cost
                _arrivalStage[arrivalState].update(_exitStage[iState],symbolVectors.col(channel.memory()-1),newCost);
    } // for(int iInput=0;iInput<_trellis->NpossibleInputs();iInput++)

}

MatrixXd ViterbiAlgorithm::getDetectedSymbolVectors()
{
    return _detectedSymbolVectors->block(0,_preamble.cols(),_nInputs,_iLastSymbolVectorToBeDetected-_preamble.cols());
}

void ViterbiAlgorithm::PrintStage(tStage exitOrArrival)
{
    ViterbiPath *stage;
    if(exitOrArrival == exitStage)
        stage = _exitStage;
    else
        stage = _arrivalStage;

    for(int i=0;i<_trellis->Nstates();i++)
    {
        cout << "State " << i << endl;
        if(stage[i].IsEmpty())
            cout << "Empty" << endl;
        else
		stage[i].print();
        cout << "------------------" << endl;
    }
}

void ViterbiAlgorithm::swapStages()
{
  // _arrivalStage becomes _exitStage for the next iteration
  ViterbiPath *aux = _exitStage;
  _exitStage = _arrivalStage;
  _arrivalStage = aux;

  // the _arrivalStage (old _exitStage) gets cleaned
  for(int iState=0;iState<_trellis->Nstates();iState++)
	  _arrivalStage[iState].clean();
}
