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

Algorithm::Algorithm(string name, Alphabet  alphabet):_name(name),_alphabet(alphabet)
{
}

double Algorithm::SER(const tMatrix &symbols)
{
    int windowSize = symbols.cols();

    tMatrix detectedSymbolVectors = GetDetectedSymbolVectors();
    int nDetectedVectors = detectedSymbolVectors.cols();

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

    // if the algorithm didn't make channel estimation
    if(nEstimatedChannelMatrices==0)
        return 0.0;

    if(windowSize>nEstimatedChannelMatrices)
        throw RuntimeException("Algorithm::MMSE: more channel matrices passed than detected.");

    double mse = 0;
    int windowStart = nEstimatedChannelMatrices - windowSize;
    int j;

    for(int i=windowStart;i<nEstimatedChannelMatrices;i++)
    {
        // the square error committed by the estimated matrix is normalized by the squared Frobenius norm (i.e. the sum of all the elements squared) of the real channel matrix
        mse += Util::SquareError(channelMatrices.at(i-windowStart),estimatedChannelMatrices.at(i))/pow(Blas_NormF(channelMatrices.at(i-windowStart)),2.0);

        cout << "sumando: " << Util::SquareError(channelMatrices.at(i-windowStart),estimatedChannelMatrices.at(i))/pow(Blas_NormF(channelMatrices.at(i-windowStart)),2.0) << ",";
    }
    cout << endl;

    return mse/(double)windowSize;
}
