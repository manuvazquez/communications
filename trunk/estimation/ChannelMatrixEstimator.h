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
#include <Util.h>

class ChannelMatrixEstimator{
protected:
    int _nOutputs,_nChannelMatrixRows,_nInputsXchannelOrder,_nInputs,_channelOrder,_nChannelCoeffs;
    MatrixXd _lastEstimatedChannelMatrix_eigen;

public:
    // initialEstimation is basically what LastEstimatedChannelMatrix is going to return when NextMatrix hasn't yet been called
    ChannelMatrixEstimator(MatrixXd initialEstimation,int N);
    virtual ~ChannelMatrixEstimator() {};

    virtual void setFirstEstimatedChannelMatrix(const MatrixXd &matrix) { _lastEstimatedChannelMatrix_eigen = matrix;} // eigen
    
    virtual MatrixXd nextMatrix(const VectorXd &observations,const MatrixXd &symbolsMatrix,double noiseVariance) = 0;
    
    virtual ChannelMatrixEstimator *clone() const = 0;
    
    virtual double likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance)
    {
        throw RuntimeException("ChannelMatrixEstimator::likelihood: not implemented yet.");
    }
    
    int cols() { return _nInputsXchannelOrder;}
    int rows() { return _nOutputs;}
    int memory();
    virtual MatrixXd lastEstimatedChannelMatrix_eigen() { return _lastEstimatedChannelMatrix_eigen;}
    
    vector<MatrixXd> nextMatricesFromObservationsSequence(const MatrixXd &observations,vector<double> &noiseVariances,const MatrixXd &symbolVectors,int iFrom,int iTo);
};

#endif
