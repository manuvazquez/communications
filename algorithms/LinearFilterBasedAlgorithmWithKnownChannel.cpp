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


#include "LinearFilterBasedAlgorithmWithKnownChannel.h"

LinearFilterBasedAlgorithmWithKnownChannel::LinearFilterBasedAlgorithmWithKnownChannel(std::string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, LinearDetector* linearDetector, std::vector< double, std::allocator< double > > ARcoefficients, bool substractContributionFromKnownSymbols): LinearFilterBasedAlgorithm(name, alphabet, L, Nr, N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, linearDetector, ARcoefficients, substractContributionFromKnownSymbols)
{
}

std::vector< MatrixXd, std::allocator< MatrixXd > > LinearFilterBasedAlgorithmWithKnownChannel::getChannelMatricesToStackForSmoothing(std::vector< MatrixXd, std::allocator< MatrixXd > > ARmatricesBuffer) const
{
	std::vector<MatrixXd> res(_d+1);
	
	ChannelMatrixEstimator *channelEstimatorClone = _channelEstimator->clone();
	
	for(uint iSmoothing=0;iSmoothing<(_d+1);iSmoothing++)
		// for the parameters of "nextMatrix" anything goes since the channel "estimator" knows the channel and hence it SHOUDN'T use them
		res[iSmoothing] = channelEstimatorClone->nextMatrix(VectorXd::Zero(1,1),MatrixXd::Zero(1,1),0.0);
	
	delete channelEstimatorClone;
	
	return res;
}
