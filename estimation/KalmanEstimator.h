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
    uint _nExtStateVectorCoeffs;
    
    virtual MatrixXd buildMeasurementMatrix(const VectorXd &symbolsVector);
public:
    KalmanEstimator(const MatrixXd &initialEstimation,const MatrixXd &variances,uint N,vector<double> ARcoefficients,double ARvariance);
    KalmanEstimator(const KalmanEstimator &kalmanEstimator);
    ~KalmanEstimator();
    
    virtual MatrixXd nextMatrix(const VectorXd &observations,const MatrixXd &symbolsMatrix,double noiseVariance);
    double likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance);
    virtual KalmanEstimator *clone() const;
    virtual MatrixXd sampleFromPredictive() const; // eigen
    virtual MatrixXd getPredictive() const {return Util::toMatrix(_kalmanFilter->predictiveMean_eigen().tail(_nChannelCoeffs),rowwise,_nChannelMatrixRows); }
    virtual void setFirstEstimatedChannelMatrix(const MatrixXd &matrix);
	
	//! it returns the corresponding covariance AS STORED by the internal Kalman Filter
    virtual MatrixXd getFilteredCovariance() const {return _kalmanFilter->filteredCovariance_eigen();}
};

#endif
