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
    uint _nOutputs,_nChannelMatrixRows,_nInputsXchannelOrder,_nInputs,_channelOrder,_nChannelCoeffs;

	// this stores the last estimated channel coefficients
    MatrixXd _lastEstimatedChannelCoefficientsMatrix;

public:
	/*!
	  It builds a \ref ChannelMatrixEstimator object
	  \param initialEstimation a matrix representing the initial estimation. It's what \ref lastEstimatedChannelMatrix returns when \ref nextMatrix hasn't been called yet
	  \param N number of inputs of the system
	*/
    ChannelMatrixEstimator(MatrixXd initialEstimation,uint N);
    virtual ~ChannelMatrixEstimator() {};

    virtual void setFirstEstimatedChannelMatrix(const MatrixXd &matrix) { _lastEstimatedChannelCoefficientsMatrix = matrix;} // eigen

    /*!
	  It updates the channel estimation of this estimator (it thus modifies the state of the estimator)
	  \param observations a vector with the new observations
	  \param symbolsMatrix the symbols involved in the observations
	  \param noiseVariance the noise variance at the relevant time instant
	  \return the updated estimation of the channel
	*/
    virtual MatrixXd nextMatrix(const VectorXd &observations,const MatrixXd &symbolsMatrix,double noiseVariance) = 0;
    
    virtual ChannelMatrixEstimator *clone() const = 0;
    
    virtual double likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance)
    {
        throw RuntimeException("ChannelMatrixEstimator::likelihood: not implemented yet.");
    }
    
    uint cols() const { return _nInputsXchannelOrder;}
    uint rows() const { return _nOutputs;}
    uint memory() const;

	/*!
	  It returns the last estimated channel matrix, that is, the one that multiplied by the symbols vector gives rise to the observations. This doesnt' necessarily coincide with the matrix ONLY containing channel coefficients (though usually, it does), which is returned by \ref lastEstimatedChannelCoefficientsMatrix
	  \return the last estimated channel matrix
	*/
    virtual MatrixXd lastEstimatedChannelMatrix() const { return lastEstimatedChannelCoefficientsMatrix();}

    virtual MatrixXd lastEstimatedChannelCoefficientsMatrix() const { return _lastEstimatedChannelCoefficientsMatrix;}
    
    vector<MatrixXd> nextMatricesFromObservationsSequence(const MatrixXd &observations,vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo);
};

#endif
