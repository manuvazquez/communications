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
#ifndef KALMANESTIMATOR_H
#define KALMANESTIMATOR_H

#include <ChannelMatrixEstimator.h>

/**
    @author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <Util.h>
#include <KalmanFilter.h>
#include <StatUtil.h>

#include <Eigen/Cholesky>

class KalmanEstimator : public ChannelMatrixEstimator
{
protected:
    KalmanFilter *_kalmanFilter;
    const uint _nExtStateVectorCoeffs;
    
    virtual MatrixXd buildMeasurementMatrix(const VectorXd &symbolsVector);
public:
    KalmanEstimator(const MatrixXd &initialEstimation,const MatrixXd &variances,uint N,vector<double> ARcoefficients,double ARvariance);
    KalmanEstimator(const KalmanEstimator &kalmanEstimator);
    ~KalmanEstimator();
	
    virtual KalmanEstimator *clone() const;
    
    virtual MatrixXd nextMatrix(const VectorXd &observations,const MatrixXd &symbolsMatrix,double noiseVariance);
    double likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance);
	
    virtual MatrixXd samplePredicted() const;
    
    virtual void setFirstEstimatedChannelMatrix(const MatrixXd &matrix);
	
	//! they return the corresponding covariance AS STORED by the internal Kalman Filter
    virtual MatrixXd getFilteredCovariance() const {return _kalmanFilter->filteredCovariance();}
    virtual MatrixXd getPredictiveCovariance() const {return _kalmanFilter->predictiveCovariance();}
    
    virtual bool computesVariances() const { return true; }
    
    /**
	 * @brief It returns a matrix the same size as the channel with each coefficient representing the variance of the corresponding estimated mean (thus it does NOT return covariances)
	 *
	 * @return MatrixXd
	 **/
	virtual MatrixXd getVariances() const 
	{
		return Util::toMatrix(_kalmanFilter->filteredCovariance().bottomRightCorner(_nChannelCoeffs,_nChannelCoeffs).diagonal(),rowwise,_nChannelMatrixRows);
	}
	
	virtual MatrixXd predictedMatrix() const { return Util::toMatrix(_kalmanFilter->predictiveMean().tail(_nChannelCoeffs),rowwise,_nChannelMatrixRows); }
	
	/**
	 * @brief It returns the indexes of the elements within the internal Kalman Filter state vector that are associated with a certain column
	 *
	 * @param iCol the index of the column
	 * @return a vector of indexes
	 **/
	std::vector<uint> colIndexToIndexesWithinKFstateVector(uint iCol) const;
};

#endif
