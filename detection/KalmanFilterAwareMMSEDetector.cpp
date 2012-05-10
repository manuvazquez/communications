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


#include "KalmanFilterAwareMMSEDetector.h"

KalmanFilterAwareMMSEDetector::KalmanFilterAwareMMSEDetector(uint rows, uint cols, double alphabetVariance,uint nSymbolsToBeDetected,KalmanEstimator *kalmanEstimator)
:MMSEDetector(rows,cols,alphabetVariance,nSymbolsToBeDetected),_kalmanEstimator(kalmanEstimator)
{
}

VectorXd KalmanFilterAwareMMSEDetector::detect(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance)
{
	uint nRows = channelMatrix.rows();
	
	MatrixXd aux = MatrixXd::Zero(nRows,nRows);
	
	for(uint i=0;i<channelMatrix.cols();i++)
		aux += covarianceMatrixForCol(i) + _kalmanEstimator->predictedMatrix().col(i)*(_kalmanEstimator->predictedMatrix().col(i)).transpose();
	
//     MatrixXd _Rx = noiseCovariance + _alphabetVariance*channelMatrix*channelMatrix.transpose();
	MatrixXd _Rx = noiseCovariance + _alphabetVariance*aux;

    _filter = _Rx.inverse()*channelMatrix*_alphabetVariance;

    VectorXd softEstimations = _filter.transpose()*observations;

    // required for nthSymbolVariance computing
    _channelMatrix = channelMatrix;

    return softEstimations.segment(_detectionStart,_nSymbolsToBeDetected);
}

KalmanFilterAwareMMSEDetector* KalmanFilterAwareMMSEDetector::clone()
{
	return new KalmanFilterAwareMMSEDetector(*this);
}

MatrixXd KalmanFilterAwareMMSEDetector::covarianceMatrixForCol(uint iCol) const
{
	MatrixXd overallCovariance = _kalmanEstimator->getPredictiveCovariance();
	
	uint n = overallCovariance.rows();
	uint nCols = _kalmanEstimator->cols();
	uint nRows = _kalmanEstimator->rows();
	
	MatrixXd res(nRows,nRows);
	
	for(uint iRow=0;iRow<nRows;iRow++)
		for(uint i=iCol,iRes=0;i<n;i+=nCols,iRes++)
			res(iRow,iRes) = overallCovariance(iRow*nCols+iCol,i);
	
	return res;
}
