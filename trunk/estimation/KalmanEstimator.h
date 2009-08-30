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
    int _nExtStateVectorCoeffs;
    
//     virtual tMatrix buildMeasurementMatrix(const tVector &symbolsVector) { return Util::eigen2lapack(buildMeasurementMatrix(Util::lapack2eigen(symbolsVector))); }
    virtual MatrixXd buildMeasurementMatrix(const VectorXd &symbolsVector); // eigen
public:
    KalmanEstimator(const MatrixXd &initialEstimation,const MatrixXd &variances,int N,vector<double> ARcoefficients,double ARvariance);
    KalmanEstimator(const KalmanEstimator &kalmanEstimator);
    ~KalmanEstimator();
    
    virtual MatrixXd nextMatrix(const VectorXd &observations,const MatrixXd &symbolsMatrix,double noiseVariance);
//     double likelihood(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance);
    double likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance);
//     {
//         return likelihood(Util::eigen2lapack(observations),Util::eigen2lapack(symbolsMatrix),noiseVariance);
//     }
    virtual KalmanEstimator *clone() const;
    virtual MatrixXd sampleFromPredictive_eigen() const; // eigen
    virtual void setFirstEstimatedChannelMatrix(const MatrixXd &matrix); // eigen
};

#endif
