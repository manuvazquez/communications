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


#include "CDMARLSEstimator.h"


CDMARLSEstimator::CDMARLSEstimator(MatrixXd initialEstimation,uint N,double forgettingFactor,const MatrixXd &spreadingCodes):RLSEstimator(initialEstimation,N,forgettingFactor),
_forgettingFactor(forgettingFactor),_R(MatrixXd::Zero(N,N)),_p(VectorXd::Zero(N)),_spreadingCodes(spreadingCodes)
{
}

MatrixXd CDMARLSEstimator::nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance)
{
	// NOTE: this code is not valid if complex symbols are allowed
	_R += _forgettingFactor*_R + symbolsMatrix.asDiagonal()*_spreadingCodes.transpose()*_spreadingCodes*symbolsMatrix.asDiagonal();
	_p += _forgettingFactor*_p + symbolsMatrix.asDiagonal()*_spreadingCodes.transpose()*observations;

	// PartialPivLU (though it's supposed to be much faster) gives rise to some NaN's
	FullPivLU<MatrixXd> luforR(_R);
	
	_lastEstimatedChannelCoefficientsMatrix = (luforR.solve(_p)).transpose();
	
	return lastEstimatedChannelMatrix();
}

MatrixXd CDMARLSEstimator::lastEstimatedChannelMatrix() const
{
	return _spreadingCodes*_lastEstimatedChannelCoefficientsMatrix.asDiagonal();
}

