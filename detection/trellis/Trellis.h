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
#ifndef TRELLIS_H
#define TRELLIS_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <Alphabet.h>
#include <math.h>

class Trellis{
protected:
    uint _nStates,_nPossibleInputs;
    uint **_stateTransitionMatrix;
public:
    Trellis(const Alphabet &alphabet, uint N, uint m);

    ~Trellis();

    uint nStates() const { return _nStates;}
    uint nPossibleInputs() const {return _nPossibleInputs;}
    uint operator ()(int state, int input) const { return _stateTransitionMatrix[state][input];}

};

#endif
