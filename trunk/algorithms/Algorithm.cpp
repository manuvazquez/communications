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

// #define DEBUG2

Algorithm::Algorithm(string name, Alphabet  alphabet,int L,int N,int K):_name(name),_alphabet(alphabet),_L(L),_N(N),_K(K)
{
}

double Algorithm::SER(const tMatrix &symbols)
{
    int windowSize = symbols.cols();

    tMatrix detectedSymbolVectors = GetDetectedSymbolVectors();
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

    vector<tMatrix> estimatedChannelMatrices = GetEstimatedChannelMatrices();
    int nEstimatedChannelMatrices = estimatedChannelMatrices.size();

    #ifdef DEBUG
    	cout << "recibidas: " << windowSize << ", estimadas: " << nEstimatedChannelMatrices << endl;
    #endif

    // if the algorithm didn't make channel estimation
    if(nEstimatedChannelMatrices==0)
        return 0.0;

    if(windowSize>nEstimatedChannelMatrices)
        throw RuntimeException("Algorithm::MMSE: more channel matrices passed than detected.");

    double mse = 0;
    int windowStart = nEstimatedChannelMatrices - windowSize;

	// if the channel is Sparkling memory, the channel matrices of the real channel may have different sizes
	try {
		for(int i=windowStart;i<nEstimatedChannelMatrices;i++)
		{
			// the square error committed by the estimated matrix is normalized by the squared Frobenius norm (i.e. the sum of all the elements squared) of the real channel matrix
			mse += Util::SquareErrorPaddingWithZeros(channelMatrices.at(i-windowStart),estimatedChannelMatrices.at(i))/pow(Blas_NormF(channelMatrices.at(i-windowStart)),2.0);
		}
	} catch (IncompatibleOperandsException) {
		return 0.0;
	}

    return mse/(double)windowSize;
}


tMatrix Algorithm::HsToStackedH(vector<tMatrix> matrices,int m,int start,int d)
{
	if((matrices[0].cols() % m)!=0)
		throw RuntimeException("Algorithm::HsToStackedH: Incorrect number of columns in the matrices.");

	uint nMatricesToStack = d - start + 1;

	if(matrices.size()< nMatricesToStack)
		throw RuntimeException("Algorithm::HsToStackedH: insufficient number of matrices.");

	tMatrix res(_L*nMatricesToStack,_N*(m+nMatricesToStack-1));
    res = 0.0;

	int iShifted;
	for(int i=start;i<=d;i++)
	{
		iShifted = i - start;
		tRange rowsRange(iShifted*_L,(iShifted+1)*_L-1);
		tRange colsRange(iShifted*_N,iShifted*_N+_N*m-1);
		res(rowsRange,colsRange).inject(matrices[i]);
	}

	return res;
}

tMatrix Algorithm::SubstractingChannelMatrix(const vector<tMatrix> &matrices,int m,int c,int d)
{
    if(matrices.size()!=c+d+1)
      throw RuntimeException("Algorithm::SubstractingChannelMatrix: wrong number of matrices.");

    int i;
    tRange rAllObservationsRows(0,_L-1);
    tMatrix res = LaGenMatDouble::zeros(_L*(c+d+1),_N*(m-1+c));

    tRange rRows(0,_L-1);
    tRange rCols(0,_N*m-1);

    for(i=0;i<c;i++)
    {
        res(rRows,rCols).inject(matrices[i]);

        rRows = rRows + _L;
        rCols = rCols + _N;
    }

#ifdef DEBUG2
    cout << "mitad" << endl;
#endif

    tRange rSourceCols(0,(m-1)*_N-1);
    int rSourceColsEnd = (m-1)*_N-1;
    for(i=c;i<c+m-1;i++)
    {
      rCols.set(_N*i,(c+m-1)*_N-1);
#ifdef DEBUG2
    cout << rRows << endl << rCols << endl << rSourceCols << endl;
#endif
      res(rRows,rCols).inject(matrices[i](rAllObservationsRows,rSourceCols));

      rRows = rRows + _L;
      rSourceColsEnd -= _N;
      rSourceCols.set(0,rSourceColsEnd);
    }

    return res;

//             for(iSmoothing=0;iSmoothing<_d;iSmoothing++)
//             {
//                 stackedChannelMatrixSubstract(tRange(iSmoothing*_L,(iSmoothing+1)*_L-1),tRange(_N*iSmoothing,(_m-1)*_N-1)).inject(matricesToStack[iSmoothing](rAllObservationsRows,tRange(0,(_m-iSmoothing-1)*_N-1)));
//             }

}
