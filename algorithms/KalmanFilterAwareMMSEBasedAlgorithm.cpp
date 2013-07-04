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


#include "KalmanFilterAwareMMSEBasedAlgorithm.h"

KalmanFilterAwareMMSEBasedAlgorithm::KalmanFilterAwareMMSEBasedAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, uint m, KalmanEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, KalmanFilterAwareMMSEDetector* linearDetector, std::vector< double, std::allocator< double > > ARcoefficients, bool substractContributionFromKnownSymbols): LinearFilterKFBasedAlgorithm(name, alphabet, L, Nr, N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, linearDetector, ARcoefficients, substractContributionFromKnownSymbols)
{
	(dynamic_cast<KalmanFilterAwareMMSEDetector *> (_linearDetector))->setKalmanEstimator(dynamic_cast<KalmanEstimator *>(_channelEstimator));
}

