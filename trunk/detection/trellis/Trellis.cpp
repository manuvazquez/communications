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
#include "Trellis.h"

Trellis::Trellis(const Alphabet &alphabet, int N, int m)
{
    _nStates = (int)pow((double)alphabet.length(),N*(m-1));
    _nPossibleInputs = (int)pow((double)alphabet.length(),N);

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

Trellis::~Trellis()
{
    for(int iState=0;iState<_nStates;iState++)
        delete[] _stateTransitionMatrix[iState];

    delete[] _stateTransitionMatrix;
}


