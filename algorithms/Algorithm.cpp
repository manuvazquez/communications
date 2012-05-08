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

MatrixXd Algorithm::channelMatrices2stackedChannelMatrix(vector< MatrixXd > matrices, uint m, uint start, uint d)
{
    if((matrices[0].cols() % m)!=0)
        throw RuntimeException("Algorithm::channelMatrices2stackedChannelMatrix: incorrect number of columns in the matrices.");

    uint nMatricesToStack = d - start + 1;

    if(matrices.size()< nMatricesToStack)
        throw RuntimeException("Algorithm::channelMatrices2stackedChannelMatrix: insufficient number of matrices.");

    MatrixXd res = MatrixXd::Zero(_nOutputs*nMatricesToStack,_nInputs*(m+nMatricesToStack-1));

    uint iStartingFromZero;
    for(uint i=start;i<=d;i++)
    {
        iStartingFromZero = i - start;
        res.block(iStartingFromZero*_nOutputs,iStartingFromZero*_nInputs,_nOutputs,_nInputs*m) = matrices[i];
    }

    return res;
}

VectorXd Algorithm::substractKnownSymbolsContribution(const vector<MatrixXd> &matrices,uint m,uint d,const VectorXd &observations,const MatrixXd &involvedSymbolVectors)
{	
	assert(matrices.size()==d+1);
	assert(observations.size()==(_nOutputs*(d+1)));
	assert(m>0);
	
    if(involvedSymbolVectors.cols()!=m-1)
         throw RuntimeException("Algorithm::substractKnownSymbolsContribution: wrong number of symbol vectors.");

    uint i;
    MatrixXd substractingChannelMatrix = MatrixXd::Zero(_nOutputs*(d+1),_nInputs*(m-1));

    for(i=0;i<m-1;i++)
        substractingChannelMatrix.block(i*_nOutputs,_nInputs*i,_nOutputs,(m-1-i)*_nInputs) = matrices[i].block(0,0,_nOutputs,(m-1-i)*_nInputs);

    return observations - substractingChannelMatrix*Util::toVector(involvedSymbolVectors,columnwise);
}
