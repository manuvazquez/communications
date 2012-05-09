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


#include "LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator.h"

LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator::LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator(const MatrixXd &initialEstimation,const MatrixXd &variances,uint N,vector<double> ARcoefficients,double ARvariance,LinearDetector *linearDetector,double symbolsVariance)
:KalmanEstimator(initialEstimation,variances,N,ARcoefficients,ARvariance),_linearDetector(linearDetector),_symbolsVariance(symbolsVariance)
{

}

LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator* LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator::clone() const
{
	return new LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator(*this);
}

MatrixXd LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator::nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance)
{
	return KalmanEstimator::nextMatrix(observations, symbolsMatrix, noiseVariance + computeExtraVariance(noiseVariance));
}

double LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator::likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance)
{
	return KalmanEstimator::likelihood(observations, symbolsMatrix, noiseVariance + computeExtraVariance(noiseVariance));
}


double LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator::computeExtraVariance(double noiseVariance)
{
	double expectationSoftEstTimesTrueSymbol = 0.0;
	double expectationSoftEstTimesSoftEst = 0.0;
	
	MatrixXd filter = _linearDetector->computedFilter();
	
	// if no filter has been computed yet (e.g., during the training sequence), the variance is not adjusted
	if(filter.cols()==0 || filter.rows()==0)
		return 0.0;
	
	expectationSoftEstTimesSoftEst = (filter.transpose()*(_symbolsVariance*lastEstimatedChannelMatrix()*lastEstimatedChannelMatrix().transpose() + noiseVariance*MatrixXd::Identity(filter.rows(),filter.rows()))*filter).trace();
	expectationSoftEstTimesTrueSymbol = (filter.transpose()*lastEstimatedChannelMatrix()).trace();
	
	double extraVariance = expectationSoftEstTimesSoftEst - 2*expectationSoftEstTimesTrueSymbol + filter.cols()*_symbolsVariance;
	
	return extraVariance;
}