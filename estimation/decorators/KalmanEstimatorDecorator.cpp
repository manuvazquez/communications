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


#include "KalmanEstimatorDecorator.h"

KalmanEstimatorDecorator::KalmanEstimatorDecorator(KalmanEstimator *kalmanEstimator)
{
	_decorated = kalmanEstimator->clone();
// 	_kalmanFilter = NULL; // this will be "deleted" by the super class!!
}

KalmanEstimatorDecorator::KalmanEstimatorDecorator(const KalmanEstimatorDecorator& other):_decorated(new KalmanEstimator(*(other._decorated)))
{
// 	_kalmanFilter = NULL; // this will be "deleted" by the super class!!
}

KalmanEstimatorDecorator::~KalmanEstimatorDecorator()
{
	delete _decorated;
}

KalmanEstimatorDecorator* KalmanEstimatorDecorator::clone() const
{
	return new KalmanEstimatorDecorator(*this);
}

vector<MatrixXd> KalmanEstimatorDecorator::nextMatricesFromObservationsSequence(const MatrixXd &observations,vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo)
{
	// since _decorated is a KF (and NOT a "decorated" KF), this will actually call KalmanEstimator's "nextMatrix" method...which is the desired action
	return _decorated->nextMatricesFromObservationsSequence(observations,noiseVariances,symbolVectors,iFrom,iTo);
}

std::vector<MatrixXd> KalmanEstimatorDecorator::nextMatricesFromObservationsSequence(const MatrixXd &observations,std::vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo,std::vector<MatrixXd> &channelEstimatesVariances)
{
	// since _decorated is a KF (and NOT a "decorated" KF), this will actually call KalmanEstimator's "nextMatrix" method...which is the desired action
	return _decorated->nextMatricesFromObservationsSequence(observations,noiseVariances,symbolVectors,iFrom,iTo,channelEstimatesVariances);
}
