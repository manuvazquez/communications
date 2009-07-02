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
#include "ChannelMatrixEstimator.h"

ChannelMatrixEstimator::ChannelMatrixEstimator(int N):_L(0),_Nm(0),_N(N),_nChannelCoeffsToBeEstimated(0)
{
}

ChannelMatrixEstimator::ChannelMatrixEstimator(tMatrix initialEstimation,int N):_L(initialEstimation.rows()),_Nm(initialEstimation.cols()),_N(N),_lastEstimatedChannelMatrix(initialEstimation),_nChannelCoeffsToBeEstimated(initialEstimation.rows()*initialEstimation.cols())
{
    if(_Nm < _N)
        throw RuntimeException("ChannelMatrixEstimator::ChannelMatrixEstimator: number of columns of \"initialEstimation\"  is less than N");

    // check erased because of "OneChannelOrderPerTransmitAtennaWrapperEstimator"
    /*
    if((_Nm % _N) != 0)
        throw RuntimeException("ChannelMatrixEstimator::ChannelMatrixEstimator: number of columns of \"initialEstimation\"  is not a multiple of N");
    */

    if((_Nm % _N) == 0)
        _m = _Nm/_N;
    // _m=-1 accounts for the case of a "OneChannelOrderPerTransmitAtennaWrapperEstimator" being used, whose internal ChannelMatrixEstimator does not need to have a number of columns multiple of _N
    else
        _m = -1;
}

int ChannelMatrixEstimator::memory()
{
    if(_m!=-1)
        return _m;
    else
        throw RuntimeException("ChannelMatrixEstimator::Memory: this may not be a real channel matrix estimator: its number of columns is not a multiple of the number of transmitting antennas.");
}

vector<tMatrix> ChannelMatrixEstimator::nextMatricesFromObservationsSequence(const tMatrix &observations,vector<double> &noiseVariances,const tMatrix &symbolVectors,int iFrom,int iTo)
{
//     tMatrix toProcessSequence = Util::append(_preamble,trainingSequence);
//     int lengthToProcessSequence = toProcessSequence.cols();

    if(observations.cols()<iTo)
        throw RuntimeException("ChannelMatrixEstimator::NextMatricesFromObservationsSequence: insufficient number of observations.");

    vector<tMatrix> estimatedMatrices(iTo-iFrom);

    // selects all the rows from a symbols matrix
    tRange allSymbolRows;

    tRange mColumns(iFrom-_m+1,iFrom);

    for(int i=iFrom;i<iTo;i++)
    {
        estimatedMatrices[i-iFrom] = nextMatrix(observations.col(i),symbolVectors(allSymbolRows,mColumns),noiseVariances[i]);
        mColumns = mColumns + 1;
    }
    return estimatedMatrices;
}
