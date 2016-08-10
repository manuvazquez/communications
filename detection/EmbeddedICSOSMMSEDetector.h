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


#ifndef EMBEDDEDICSOSMMSEDETECTOR_H
#define EMBEDDEDICSOSMMSEDETECTOR_H

#include <SOSMMSEDetector.h>

#include<iostream>
#include <KalmanEstimator.h>

class EmbeddedICSOSMMSEDetector : public SOSMMSEDetector
{
	
public:
    EmbeddedICSOSMMSEDetector(uint rows, uint cols, double alphabetVariance,uint nSymbolsToBeDetected,KalmanEstimator *kalmanEstimator,std::vector<double> ARcoefficients);
	
    virtual VectorXd detect(const VectorXd &observations, const MatrixXd &channelMatrix, const MatrixXd& noiseCovariance);
	
    virtual EmbeddedICSOSMMSEDetector* clone();
};

#endif // EMBEDDEDICSOSMMSEDETECTOR_H
