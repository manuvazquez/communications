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


#include "LinearFilterKFBasedAlgorithm.h"

std::vector<MatrixXd> LinearFilterKFBasedAlgorithm::getChannelMatricesToStackForSmoothing(std::vector<MatrixXd> ARmatricesBuffer) const
{
	std::vector<MatrixXd> res(_d+1);
	res[0] = _channelEstimator->predictedMatrix();
	
	KalmanEstimator *kalmanEstimatorClone = dynamic_cast<KalmanEstimator *>(_channelEstimator)->clone();
	
	for(uint i=1;i<=_d;i++)
	{
		// one PREDICTIVE step is advanced on the cloned Kalman filter: observations vector and symbols matrix are zero so that the filtered mean and covariance be equal to the predictive ones.
		// As for the noise variance, a value different from 0 (any should do) must be passed or otherwise, the matrix inversion within the KF will give rise to NaN's
		kalmanEstimatorClone->nextMatrix(VectorXd::Zero(_nOutputs),MatrixXd::Zero(_nInputs,_channelOrder),1.0);
		
		res[i] = kalmanEstimatorClone->predictedMatrix();
	}
	
	delete kalmanEstimatorClone;
	
	return res;
}

LinearFilterKFBasedAlgorithm::LinearFilterKFBasedAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, KalmanEstimator* kalmanEstimator, MatrixXd preamble, uint smoothingLag, LinearDetector *linearDetector, std::vector<double> ARcoefficients, bool substractContributionFromKnownSymbols): LinearFilterBasedAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, kalmanEstimator, preamble, smoothingLag, linearDetector, ARcoefficients, substractContributionFromKnownSymbols)
{

}

