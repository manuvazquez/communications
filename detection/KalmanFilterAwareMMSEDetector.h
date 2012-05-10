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


#ifndef KALMANFILTERAWAREMMSEDETECTOR_H
#define KALMANFILTERAWAREMMSEDETECTOR_H

#include <MMSEDetector.h>


#include <KalmanEstimator.h>

class KalmanFilterAwareMMSEDetector : public MMSEDetector
{
protected:
	KalmanEstimator const *_kalmanEstimator;
public:
    KalmanFilterAwareMMSEDetector(uint rows, uint cols, double alphabetVariance,uint nSymbolsToBeDetected,KalmanEstimator *kalmanEstimator);
	
    virtual VectorXd detect(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance);
	
    virtual KalmanFilterAwareMMSEDetector* clone();
	
	MatrixXd covarianceMatrixForCol(uint iCol) const;
	
	void setKalmanEstimator(KalmanEstimator * kalmanEstimator) { _kalmanEstimator = kalmanEstimator; }
};

#endif // KALMANFILTERAWAREMMSEDETECTOR_H
