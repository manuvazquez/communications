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

// #define DEBUG

Algorithm::Algorithm(string name, Alphabet  alphabet,int L,int Nr,int N,int iLastSymbolVectorToBeDetected):_name(name),_alphabet(alphabet),_nOutputs(L),_Nr(Nr),_nInputs(N),_iLastSymbolVectorToBeDetected(iLastSymbolVectorToBeDetected)
{
}

double Algorithm::SER(const tMatrix &symbols)
{
    int windowSize = symbols.cols();

    tMatrix detectedSymbolVectors = getDetectedSymbolVectors();
    int nDetectedVectors = detectedSymbolVectors.cols();

    // if the algorithm did know the symbols
    if(detectedSymbolVectors.cols()==0)
        return 0.0;

    if(windowSize>nDetectedVectors)
        throw RuntimeException("Algorithm::SER: more symbol vectors passed than detected.");

    int nErrors = 0;
    int windowStart = nDetectedVectors - windowSize;
    int j;

    for(int i=windowStart;i<nDetectedVectors;i++)
    {
        j=0;
        while(j<symbols.rows())
        {
            if(detectedSymbolVectors(j,i)!=symbols(j,i-windowStart))
                nErrors++;
            j++;
        }
    }
    return ((double)nErrors)/(double)(windowSize*symbols.rows());
}

double Algorithm::MSE(const vector<tMatrix> &channelMatrices)
{
    int windowSize = channelMatrices.size();

    vector<tMatrix> estimatedChannelMatrices = getEstimatedChannelMatrices();
    int nEstimatedChannelMatrices = estimatedChannelMatrices.size();

    #ifdef DEBUG
        cout << "recibidas: " << windowSize << ", estimadas: " << nEstimatedChannelMatrices << endl;
    #endif

    // if the algorithm didn't make channel estimation
    if(nEstimatedChannelMatrices==0)
        return 0.0;

    if(windowSize>nEstimatedChannelMatrices)
        throw RuntimeException("Algorithm::MSE: more channel matrices passed than detected.");

    double mse = 0;
    int windowStart = nEstimatedChannelMatrices - windowSize;

    // if the channel is Sparkling memory, the channel matrices of the real channel may have different sizes
    try {
        for(int i=windowStart;i<nEstimatedChannelMatrices;i++)
        {
#ifdef DEBUG
            cout << "channelMatrices.at(i-windowStart) = " << endl << channelMatrices.at(i-windowStart);
            cout << "estimatedChannelMatrices.at(i) = " << endl << estimatedChannelMatrices.at(i);
#endif
            // the square error committed by the estimated matrix is normalized by the squared Frobenius norm (i.e. the sum of all the elements squared) of the real channel matrix
            mse += Util::squareErrorPaddingWithZeros(channelMatrices.at(i-windowStart),estimatedChannelMatrices.at(i))/pow(Blas_NormF(channelMatrices.at(i-windowStart)),2.0);
        }
    } catch (IncompatibleOperandsException) {
        return 0.0;
    }

    return mse/(double)windowSize;
}

MatrixXd Algorithm::channelMatrices2stackedChannelMatrix(vector<MatrixXd> matrices,int m,int start,int d)
{
    if((matrices[0].cols() % m)!=0)
        throw RuntimeException("Algorithm::channelMatrices2stackedChannelMatrix: incorrect number of columns in the matrices.");

    uint nMatricesToStack = d - start + 1;

    if(matrices.size()< nMatricesToStack)
        throw RuntimeException("Algorithm::channelMatrices2stackedChannelMatrix: insufficient number of matrices.");

    MatrixXd res = MatrixXd::Zero(_nOutputs*nMatricesToStack,_nInputs*(m+nMatricesToStack-1));

    int iStartingFromZero;
    for(int i=start;i<=d;i++)
    {
        iStartingFromZero = i - start;
        res.block(iStartingFromZero*_nOutputs,iStartingFromZero*_nInputs,_nOutputs,_nInputs*m) = matrices[i];
    }

    return res;
}

VectorXd Algorithm::substractKnownSymbolsContribution(const vector<MatrixXd> &matrices,int m,int c,int e,const VectorXd &observations,const MatrixXd &involvedSymbolVectors)
{
    if(matrices.size()!=static_cast<uint> (c+e+1))
      throw RuntimeException("Algorithm::substractKnownSymbolsContribution: wrong number of matrices.");

    if(observations.size()!=(_nOutputs*(c+e+1)))
       throw RuntimeException("Algorithm::substractKnownSymbolsContribution: size of observations vector is wrong.");

    if(involvedSymbolVectors.cols()!=c+m-1)
         throw RuntimeException("Algorithm::substractKnownSymbolsContribution: wrong number of symbol vectors.");

    int i;
    MatrixXd substractingChannelMatrix = MatrixXd::Zero(_nOutputs*(c+e+1),_nInputs*(m-1+c));

    for(i=0;i<c;i++)
        substractingChannelMatrix.block(i*_nOutputs,i*_nInputs,_nOutputs,_nInputs*m) = matrices[i];

    for(i=c;i<c+m-1;i++)
        substractingChannelMatrix.block(i*_nOutputs,_nInputs*i,_nOutputs,(c+m-1-i)*_nInputs) = matrices[i].block(0,0,_nOutputs,(m-1-(i-c))*_nInputs);

    return observations - substractingChannelMatrix*Util::toVector(involvedSymbolVectors,columnwise);
}
