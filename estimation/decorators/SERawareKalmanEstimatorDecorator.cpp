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


#include "SERawareKalmanEstimatorDecorator.h"

SERawareKalmanEstimatorDecorator::SERawareKalmanEstimatorDecorator(KalmanEstimator *kalmanEstimator, double symbolsVariance):
KalmanEstimatorDecorator(kalmanEstimator),_symbolsVariance(symbolsVariance)
{
}

MatrixXd SERawareKalmanEstimatorDecorator::nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance)
{
	return KalmanEstimatorDecorator::nextMatrix(observations, symbolsMatrix, noiseVariance + computeExtraVariance(noiseVariance));
}

double SERawareKalmanEstimatorDecorator::likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance)
{
	return KalmanEstimatorDecorator::likelihood(observations, symbolsMatrix, noiseVariance + computeExtraVariance(noiseVariance));
}

SERawareKalmanEstimatorDecorator* SERawareKalmanEstimatorDecorator::clone() const
{
	return new SERawareKalmanEstimatorDecorator(*this);
}

double SERawareKalmanEstimatorDecorator::computeExtraVariance(double noiseVariance)
{
	return 0.0;
}


