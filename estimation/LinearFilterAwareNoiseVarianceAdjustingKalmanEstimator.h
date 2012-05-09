/*
    Copyright 2012 Manu <manuavazquez@gmail.com>

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/


#ifndef LINEARFILTERAWARENOISEVARIANCEADJUSTINGKALMANESTIMATOR_H
#define LINEARFILTERAWARENOISEVARIANCEADJUSTINGKALMANESTIMATOR_H

#include <KalmanEstimator.h>

#include <LinearDetector.h>

class LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator : public KalmanEstimator
{
protected:
	LinearDetector const *_linearDetector;
	double _symbolsVariance;
	
	double computeExtraVariance(double noiseVariance);
	
public:
    LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator(const MatrixXd &initialEstimation,const MatrixXd &variances,uint N,vector<double> ARcoefficients,double ARvariance,LinearDetector *linearDetector,double symbolsVariance);
	
    virtual LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator* clone() const;
	
    virtual MatrixXd nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance);
    virtual double likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance);
	
	void setLinearDetector(LinearDetector * linearDetector) { _linearDetector = linearDetector; }
};

#endif // LINEARFILTERAWARENOISEVARIANCEADJUSTINGKALMANESTIMATOR_H
