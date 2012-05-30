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


#ifndef LINEARFILTERBASEDALGORITHMWITHKNOWNCHANNEL_H
#define LINEARFILTERBASEDALGORITHMWITHKNOWNCHANNEL_H

#include <LinearFilterBasedAlgorithm.h>


class LinearFilterBasedAlgorithmWithKnownChannel : public LinearFilterBasedAlgorithm
{

protected:
    virtual std::vector< MatrixXd, std::allocator< MatrixXd > > getChannelMatricesToStackForSmoothing(std::vector< MatrixXd, std::allocator< MatrixXd > > ARmatricesBuffer) const;

public:
    LinearFilterBasedAlgorithmWithKnownChannel(std::string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, LinearDetector* linearDetector, std::vector< double, std::allocator< double > > ARcoefficients, bool substractContributionFromKnownSymbols = false);
};

#endif // LINEARFILTERBASEDALGORITHMWITHKNOWNCHANNEL_H
