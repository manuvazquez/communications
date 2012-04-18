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
#include <defines.h>
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

    virtual void setFirstEstimatedChannelMatrix(const MatrixXd &matrix) { _lastEstimatedChannelCoefficientsMatrix = matrix;}

    /*!
	  It updates the channel estimation of this estimator (it thus modifies the state of the estimator)
	  \param observations a vector with the new observations
	  \param symbolsMatrix the symbols involved in the observations
	  \param noiseVariance the noise variance at the relevant time instant
	  \return the updated estimation of the channel
	*/
    virtual MatrixXd nextMatrix(const VectorXd &observations,const MatrixXd &symbolsMatrix,double noiseVariance) = 0;
	
	/**
	 * @brief It returns an estimate of the channel at the current time instant as predicted by the estimator, i.e., prior to using the corresponding observations and symbols to update it (as opposed to "nextMatrix"). It defaults to the last estimated channel matrix
	 * @return MatrixXd
	 **/
	virtual MatrixXd predictedMatrix() const { return lastEstimatedChannelMatrix();}
	
	/**
	 * @brief It returns (if possible) a sample of the matrix returned by predictedMatrix. It defaults to simply calling the latter.
	 *
	 * @return MatrixXd
	 **/
	virtual MatrixXd samplePredicted() const { return predictedMatrix(); }
    
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
	
	/**
	 * @brief It runs the channel matrix estimator algorithm for several steps in a row
	 *
	 * @param observations a sequence of observations (possibly more than will be used)
	 * @param noiseVariances the noise variances corresponding to the observations
	 * @param symbolVectors sequence of symbol vectors
	 * @param iFrom first time instant to be accounted for
	 * @param iTo last time instant to be accounted for
	 * @param channelEstimatesVariances return parameter
	 * @return channel matrices obtained from all the channel matrix estimator steps
	 **/
	std::vector<MatrixXd> nextMatricesFromObservationsSequence(const MatrixXd &observations,std::vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo,std::vector<MatrixXd> &channelEstimatesVariances);
	
	virtual bool computesVariances() const { return false; }
	
// 	virtual MatrixXd getVariances() const {return MatrixXd::Constant(_nOutputs,_nInputsXchannelOrder,FUNNY_VALUE); }
	virtual MatrixXd getVariances() const { throw RuntimeException("ChannelMatrixEstimator::getVariances: not implemented for this estimator."); }
};

#endif
