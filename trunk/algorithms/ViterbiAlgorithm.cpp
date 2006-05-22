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

ViterbiAlgorithm::ViterbiAlgorithm(string name, Alphabet alphabet, const MIMOChannel& channel,int detectionLag): KnownChannelAlgorithm(name, alphabet, channel),_detectionLag(detectionLag)
{
	_nStates = (int)pow((double)alphabet.Length(),(double)channel.Nt()*(channel.Memory()-1));
	_nPossibleInputs = (int)pow((double)alphabet.Length(),(double)channel.Nt());
	BuildStateTransitionMatrix();
}


ViterbiAlgorithm::~ViterbiAlgorithm()
{
	for(int iState=0;iState<_nStates;iState++)
		delete[] _stateTransitionMatrix[iState];

	delete[] _stateTransitionMatrix;
}

void ViterbiAlgorithm::Run(const tMatrix &observations,vector<double> noiseVariances)
{
}

void ViterbiAlgorithm::BuildStateTransitionMatrix()
{
	_stateTransitionMatrix = new int*[_nStates];

	int alphabetLengthToTheNmMinus2 = _nStates/_nPossibleInputs;

// 	cout << "nstates = " << _nStates << "nPossibleInputs = " << _nPossibleInputs << " lo otro " << alphabetLengthToTheNmMinus2 << endl;

	int iInput;
	for(int iState=0;iState<_nStates;iState++)
	{
		_stateTransitionMatrix[iState] = new int[_nPossibleInputs];
		for(iInput=0;iInput<_nPossibleInputs;iInput++)
			// computes de next state give the current one and the input (both in decimal)
			_stateTransitionMatrix[iState][iInput] = (iState % alphabetLengthToTheNmMinus2)*_nPossibleInputs + iInput;
	}

	for(int iState=0;iState<_nStates;iState++)
	{
		for(iInput=0;iInput<_nPossibleInputs;iInput++)
			cout << _stateTransitionMatrix[iState][iInput] << " ";
		cout << endl;
	}
	
}
