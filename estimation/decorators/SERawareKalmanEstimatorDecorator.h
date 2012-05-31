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


#ifndef SERAWAREKALMANESTIMATORDECORATOR_H
#define SERAWAREKALMANESTIMATORDECORATOR_H

#include <KalmanEstimatorDecorator.h>

#include <LinearDetector.h>

class SERawareKalmanEstimatorDecorator : public KalmanEstimatorDecorator
{
protected:
	const std::vector<double> _noiseVariances;
	const std::vector<double> _SERs;
	const std::vector<tSymbol> _possibleErrors;
	
	double computeExtraVariance(double noiseVariance);
public:
	SERawareKalmanEstimatorDecorator(KalmanEstimator *kalmanEstimator,const std::vector<double> &noiseVariances,const std::vector<double> &SERs,const std::vector<double> &possibleErrors);
	
    virtual MatrixXd nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance);
    virtual double likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance);
    virtual SERawareKalmanEstimatorDecorator* clone() const;
};

#endif // SERAWAREKALMANESTIMATORDECORATOR_H
