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

// #define DEBUG

Algorithm::Algorithm(string name, Alphabet  alphabet,uint L,uint Nr,uint N,uint iLastSymbolVectorToBeDetected):_name(name),_alphabet(alphabet),_nOutputs(L),_Nr(Nr),_nInputs(N),_iLastSymbolVectorToBeDetected(iLastSymbolVectorToBeDetected)
{
}

double Algorithm::MSE(const vector<MatrixXd> &channelMatrices)
{
    uint windowSize = channelMatrices.size();

    vector<MatrixXd> estimatedChannelMatrices = getEstimatedChannelMatrices();
    uint nEstimatedChannelMatrices = estimatedChannelMatrices.size();

    // if the algorithm didn't make channel estimation
    if(nEstimatedChannelMatrices==0)
        return -1.0;

    if(windowSize>nEstimatedChannelMatrices)
        throw RuntimeException("Algorithm::MSE: more channel matrices passed than detected.");

    double mse = 0;
	
	assert(nEstimatedChannelMatrices>=windowSize);
    uint windowStart = nEstimatedChannelMatrices - windowSize;

#ifdef DEBUG
	cout << "computing MSE..." << endl;
	cout << "windowStart = " << windowStart << " nEstimatedChannelMatrices = " << nEstimatedChannelMatrices << endl;
#endif
	
    // if the channel is Sparkling memory, the channel matrices of the real channel may have different sizes
	for(uint i=windowStart;i<nEstimatedChannelMatrices;i++)
	{
		// the square error committed by the estimated matrix is normalized by the squared Frobenius norm
		// (i.e. the sum of all the elements squared) of the real channel matrix
#ifdef DEBUG
		cout << "comparing" << endl << channelMatrices.at(i-windowStart) << endl << "and" << endl << estimatedChannelMatrices.at(i) << endl;
		cout << "result = " << Util::squareErrorPaddingWithZeros(channelMatrices.at(i-windowStart),estimatedChannelMatrices.at(i))/channelMatrices.at(i-windowStart).squaredNorm() << endl;
#endif
		mse += Util::squareErrorPaddingWithZeros(channelMatrices.at(i-windowStart),estimatedChannelMatrices.at(i))/channelMatrices.at(i-windowStart).squaredNorm();
	}

    return mse/(double)windowSize;
}

double Algorithm::MSE(const vector<MatrixXd> &channelMatrices,const vector<uint> &bestPermutation,const vector<int> &bestPermutationSigns)
{
  vector<uint> realChannelMatricesPermutation = Util::computeInversePermutation(bestPermutation);

  // signs permutation is given  WITH RESPECT TO THE ESTIMATED CHANNEL MATRICES. 
  // We have to permute them according to the best permutation with respect to the real channel matrices
  vector<int> realChannelMatricesSignsPermutation = Util::applyPermutation(Util::applyPermutation(bestPermutationSigns,realChannelMatricesPermutation),realChannelMatricesPermutation);

  vector<MatrixXd> permutedChannelMatrices(channelMatrices.size());
  
  for(uint i=0;i<channelMatrices.size();i++)
	permutedChannelMatrices[i] = Util::applyPermutationOnColumns(channelMatrices[i],realChannelMatricesPermutation,realChannelMatricesSignsPermutation);

  return MSE(permutedChannelMatrices);
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

VectorXd Algorithm::substractKnownSymbolsContribution(const vector<MatrixXd> &matrices,uint m,uint c,uint e,const VectorXd &observations,const MatrixXd &involvedSymbolVectors)
{
    if(matrices.size()!=c+e+1)
      throw RuntimeException("Algorithm::substractKnownSymbolsContribution: wrong number of matrices.");

    if(observations.size()!=(_nOutputs*(c+e+1)))
       throw RuntimeException("Algorithm::substractKnownSymbolsContribution: size of observations vector is wrong.");

	assert(c+m>0);
    if(involvedSymbolVectors.cols()!=c+m-1)
         throw RuntimeException("Algorithm::substractKnownSymbolsContribution: wrong number of symbol vectors.");

    uint i;
    MatrixXd substractingChannelMatrix = MatrixXd::Zero(_nOutputs*(c+e+1),_nInputs*(m-1+c));

    for(i=0;i<c;i++)
        substractingChannelMatrix.block(i*_nOutputs,i*_nInputs,_nOutputs,_nInputs*m) = matrices[i];

    for(i=c;i<c+m-1;i++)
        substractingChannelMatrix.block(i*_nOutputs,_nInputs*i,_nOutputs,(c+m-1-i)*_nInputs) = matrices[i].block(0,0,_nOutputs,(m-1-(i-c))*_nInputs);

    return observations - substractingChannelMatrix*Util::toVector(involvedSymbolVectors,columnwise);
}
