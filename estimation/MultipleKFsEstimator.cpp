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


#include "MultipleKFsEstimator.h"

MultipleKFsEstimator::MultipleKFsEstimator(MatrixXd initialEstimation, const MatrixXd& variances, uint N, std::vector< double > ARcoefficients, double ARvariance):ChannelMatrixEstimator(initialEstimation,N),_kalmanEstimators(initialEstimation.rows())
{
	for(uint i=0;i<_kalmanEstimators.size();i++)
		_kalmanEstimators[i] = new KalmanEstimator(initialEstimation.row(i),variances.row(i),N,ARcoefficients,ARvariance);
}

MultipleKFsEstimator::MultipleKFsEstimator(const MultipleKFsEstimator& multipleKFsEstimator):ChannelMatrixEstimator(multipleKFsEstimator),_kalmanEstimators(multipleKFsEstimator._kalmanEstimators.size())
{
	for(uint i=0;i<multipleKFsEstimator._kalmanEstimators.size();i++)
		_kalmanEstimators[i] = multipleKFsEstimator._kalmanEstimators[i]->clone();
}

MultipleKFsEstimator::~MultipleKFsEstimator()
{
	for(uint i=0;i<_kalmanEstimators.size();i++)
		delete _kalmanEstimators[i];
}

MultipleKFsEstimator* MultipleKFsEstimator::clone() const
{
	return new MultipleKFsEstimator(*this);
}

MatrixXd MultipleKFsEstimator::nextMatrix(const VectorXd &observations,const MatrixXd &symbolsMatrix,double noiseVariance)
{
	MatrixXd res(_nOutputs,_nInputsXchannelOrder);
	
	for(uint i=0;i<_kalmanEstimators.size();i++)
		res.row(i) = _kalmanEstimators[i]->nextMatrix(observations.segment(i,1),symbolsMatrix,noiseVariance);
	
	return res;
}

