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


#ifndef KALMANESTIMATORDECORATOR_H
#define KALMANESTIMATORDECORATOR_H

#include <KalmanEstimator.h>


class KalmanEstimatorDecorator : public KalmanEstimator
{

protected:
	KalmanEstimator *_decorated;
	
// 	virtual MatrixXd buildMeasurementMatrix(const VectorXd& symbolsVector);


public:
    KalmanEstimatorDecorator(KalmanEstimator *kalmanEstimator);
    KalmanEstimatorDecorator(const KalmanEstimatorDecorator& other);
    virtual ~KalmanEstimatorDecorator();
	
    virtual MatrixXd nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance);
    virtual double likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance);
    virtual KalmanEstimatorDecorator* clone() const;
    virtual MatrixXd samplePredicted() const;
    virtual void setFirstEstimatedChannelMatrix(const MatrixXd& matrix);
    virtual MatrixXd getFilteredCovariance() const;
    virtual bool computesVariances() const;
    virtual MatrixXd getVariances() const;
    virtual MatrixXd predictedMatrix() const;
    virtual MatrixXd lastEstimatedChannelMatrix() const;
    virtual MatrixXd lastEstimatedChannelCoefficientsMatrix() const;
	
	virtual uint cols() const { return _decorated->cols();}
	virtual uint rows() const { return _decorated->rows();}
	virtual uint memory() const {return _decorated->memory();}
	virtual vector<MatrixXd> nextMatricesFromObservationsSequence(const MatrixXd &observations,vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo);
	virtual std::vector<MatrixXd> nextMatricesFromObservationsSequence(const MatrixXd &observations,std::vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo,std::vector<MatrixXd> &channelEstimatesVariances);
};

#endif // KALMANESTIMATORDECORATOR_H
