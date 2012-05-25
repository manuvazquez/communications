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
#include "Algorithm.h"
#include <assert.h>

Algorithm::Algorithm(std::string name, Alphabet  alphabet,uint L,uint Nr,uint N,uint iLastSymbolVectorToBeDetected):_name(name),_alphabet(alphabet),_nOutputs(L),_Nr(Nr),_nInputs(N),_iLastSymbolVectorToBeDetected(iLastSymbolVectorToBeDetected)
{
}

VectorXd Algorithm::substractKnownSymbolsContribution(const std::vector<MatrixXd> &matrices,uint m,uint d,const VectorXd &observations,const MatrixXd &involvedSymbolVectors)
{	
	assert(matrices.size()==d+1);
	assert(observations.size()==(_nOutputs*(d+1)));
	assert(m>0);

	// wrong number of symbol vectors
	assert(involvedSymbolVectors.cols()==m-1);

    uint i;
    MatrixXd substractingChannelMatrix = MatrixXd::Zero(_nOutputs*(d+1),_nInputs*(m-1));

    for(i=0;i<m-1;i++)
        substractingChannelMatrix.block(i*_nOutputs,_nInputs*i,_nOutputs,(m-1-i)*_nInputs) = matrices[i].block(0,0,_nOutputs,(m-1-i)*_nInputs);

    return observations - substractingChannelMatrix*Util::toVector(involvedSymbolVectors,columnwise);
}
