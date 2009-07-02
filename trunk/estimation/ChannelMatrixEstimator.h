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
#ifndef CHANNELMATRIXESTIMATOR_H
#define CHANNELMATRIXESTIMATOR_H

/**
    @author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <exceptions.h>
#include <vector>

class ChannelMatrixEstimator{
protected:
    int _L,_Nm,_N,_m,_nChannelCoeffsToBeEstimated;
    tMatrix _lastEstimatedChannelMatrix;

    ChannelMatrixEstimator(int N);
public:
    // initialEstimation is basically what LastEstimatedChannelMatrix is going to return when NextMatrix hasn't yet been called
    ChannelMatrixEstimator(tMatrix initialEstimation,int N);
    virtual ~ChannelMatrixEstimator() {};

    virtual void setFirstEstimatedChannelMatrix(const tMatrix &matrix) { _lastEstimatedChannelMatrix = matrix;}
    virtual tMatrix nextMatrix(const tVector &observations,const tMatrix &symbolsMatrix,double noiseVariance) = 0;
    virtual ChannelMatrixEstimator *Clone() const = 0;

    virtual double likelihood(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance)
    {
        throw RuntimeException("ChannelMatrixEstimator::Likelihood: not implemented yet.");
    }
    int cols() { return _Nm;}
    int rows() { return _L;}
    int memory();
    virtual tMatrix lastEstimatedChannelMatrix() { return _lastEstimatedChannelMatrix;}
    vector<tMatrix> nextMatricesFromObservationsSequence(const tMatrix &observations,vector<double> &noiseVariances,const tMatrix &symbolVectors,int iFrom,int iTo);
};

#endif
