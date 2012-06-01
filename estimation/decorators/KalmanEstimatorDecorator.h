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
	
public:
    KalmanEstimatorDecorator(KalmanEstimator *kalmanEstimator):_decorated(kalmanEstimator->clone()) {}
    KalmanEstimatorDecorator(const KalmanEstimatorDecorator& other):_decorated(new KalmanEstimator(*(other._decorated))) {}
    virtual ~KalmanEstimatorDecorator() { delete _decorated;}
    
    virtual MatrixXd nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance) { return _decorated->nextMatrix(observations, symbolsMatrix, noiseVariance);}
    virtual double likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance) { return _decorated->likelihood(observations, symbolsMatrix, noiseVariance);}
    virtual KalmanEstimatorDecorator* clone() const { return new KalmanEstimatorDecorator(*this);}
    virtual MatrixXd samplePredicted() const { return _decorated->samplePredicted();}
    virtual void setFirstEstimatedChannelMatrix(const MatrixXd& matrix) { _decorated->setFirstEstimatedChannelMatrix(matrix);}
    virtual MatrixXd getFilteredCovariance() const { return _decorated->getFilteredCovariance();}
	virtual MatrixXd getPredictiveCovariance() const {return _decorated->getPredictiveCovariance();}
	virtual VectorXd getPredictiveMean() const {return _decorated->getPredictiveMean();}
    virtual bool computesVariances() const { return _decorated->computesVariances();}
    virtual MatrixXd getVariances() const { return _decorated->getVariances();}
    virtual MatrixXd predictedMatrix() const { return _decorated->predictedMatrix();}
    virtual MatrixXd lastEstimatedChannelMatrix() const { return _decorated->lastEstimatedChannelMatrix();}
    virtual MatrixXd lastEstimatedChannelCoefficientsMatrix() const { return _decorated->lastEstimatedChannelCoefficientsMatrix();}
	
	virtual uint cols() const { return _decorated->cols();}
	virtual uint rows() const { return _decorated->rows();}
	virtual uint nInputs() const { return _decorated->nInputs();}
	virtual uint memory() const {return _decorated->memory();}
	
	virtual vector<MatrixXd> nextMatricesFromObservationsSequence(const MatrixXd &observations,vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo)
	{
		// since _decorated is a KF (and NOT a "decorated" KF), this will actually call KalmanEstimator's "nextMatrix" method...which is the desired action
		return _decorated->nextMatricesFromObservationsSequence(observations,noiseVariances,symbolVectors,iFrom,iTo);
	}
	
	virtual std::vector<MatrixXd> nextMatricesFromObservationsSequence(const MatrixXd &observations,std::vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo,std::vector<MatrixXd> &channelEstimatesVariances)
	{
		// since _decorated is a KF (and NOT a "decorated" KF), this will actually call KalmanEstimator's "nextMatrix" method...which is the desired action
		return _decorated->nextMatricesFromObservationsSequence(observations,noiseVariances,symbolVectors,iFrom,iTo,channelEstimatesVariances);
	}
	
	virtual std::vector<uint> colIndexToIndexesWithinKFstateVector(uint iCol) const { return _decorated->colIndexToIndexesWithinKFstateVector(iCol);}
};

#endif // KALMANESTIMATORDECORATOR_H
