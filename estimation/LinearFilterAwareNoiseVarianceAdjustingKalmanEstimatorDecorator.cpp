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


#include "LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator.h"

LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator::LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator(KalmanEstimator *kalmanEstimator, LinearDetector *linearDetector, double symbolsVariance):
KalmanEstimatorDecorator(kalmanEstimator),_linearDetector(linearDetector),_symbolsVariance(symbolsVariance)
{
}

MatrixXd LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator::nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance)
{
	return KalmanEstimatorDecorator::nextMatrix(observations, symbolsMatrix, noiseVariance + computeExtraVariance(noiseVariance));
}

double LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator::likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance)
{
	return KalmanEstimatorDecorator::likelihood(observations, symbolsMatrix, noiseVariance + computeExtraVariance(noiseVariance));
}

LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator* LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator::clone() const
{
	return new LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator(*this);
}

double LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator::computeExtraVariance(double noiseVariance)
{
	double expectationSoftEstTimesTrueSymbol = 0.0;
	double expectationSoftEstTimesSoftEst = 0.0;
	
	MatrixXd filter = _linearDetector->computedFilter();
	
	expectationSoftEstTimesSoftEst = (filter.transpose()*(_symbolsVariance*lastEstimatedChannelMatrix()*lastEstimatedChannelMatrix().transpose() + noiseVariance*MatrixXd::Identity(filter.rows(),filter.rows()))*filter).trace();
	expectationSoftEstTimesTrueSymbol = (filter.transpose()*lastEstimatedChannelMatrix()).trace();
	
	double extraVariance = expectationSoftEstTimesSoftEst - 2*expectationSoftEstTimesTrueSymbol + filter.cols()*_symbolsVariance;
	
// 	expectationSoftEstTimesSoftEst = (filter.transpose()*filter).trace()*(filter.rows()*_symbolsVariance*1 + noiseVariance);
// 	double extraVariance = expectationSoftEstTimesSoftEst + filter.rows()*_symbolsVariance;
	
	return extraVariance;
}


