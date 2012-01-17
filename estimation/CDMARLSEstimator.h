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


#ifndef CDMARLSESTIMATOR_H
#define CDMARLSESTIMATOR_H

#include <RLSEstimator.h>

#include <Eigen/LU> 
#include <Eigen/Cholesky>

class CDMARLSEstimator : public RLSEstimator
{
protected:
	double _forgettingFactor;
    MatrixXd _R;
	VectorXd _p;
	MatrixXd _spreadingCodes;
public:
    CDMARLSEstimator(MatrixXd initialEstimation,uint N,double forgettingFactor,const MatrixXd &spreadingCodes);
	
    virtual CDMARLSEstimator* clone() const { return new CDMARLSEstimator(*this); }
    virtual MatrixXd nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance);
    virtual MatrixXd lastEstimatedChannelMatrix() const;
};

#endif // CDMARLSESTIMATOR_H
