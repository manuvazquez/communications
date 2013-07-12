/*
    Copyright 2013 Manu <manuavazquez@gmail.com>

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


#ifndef MULTIPLEKFSESTIMATOR_H
#define MULTIPLEKFSESTIMATOR_H

#include <ChannelMatrixEstimator.h>

#include <KalmanEstimator.h>

class MultipleKFsEstimator : public ChannelMatrixEstimator
{
protected:
	std::vector<KalmanEstimator *> _kalmanEstimators;
public:
    MultipleKFsEstimator(MatrixXd initialEstimation, const MatrixXd &variances, uint N, std::vector<double> ARcoefficients,double ARvariance);
    MultipleKFsEstimator(const MultipleKFsEstimator& multipleKFsEstimator);
    virtual ~MultipleKFsEstimator();
	
    virtual MultipleKFsEstimator* clone() const;
    virtual MatrixXd nextMatrix(const VectorXd &observations,const MatrixXd &symbolsMatrix,double noiseVariance);
};

#endif // MULTIPLEKFSESTIMATOR_H
