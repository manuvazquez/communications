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


#ifndef LINEARFILTERKFBASEDALGORITHM_H
#define LINEARFILTERKFBASEDALGORITHM_H

#include <LinearFilterBasedAlgorithm.h>

#include <KalmanEstimator.h>

class LinearFilterKFBasedAlgorithm : public LinearFilterBasedAlgorithm
{

protected:
    virtual std::vector<MatrixXd> getChannelMatricesToStackForSmoothing(std::vector<MatrixXd> ARmatricesBuffer) const;

public:
    LinearFilterKFBasedAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, KalmanEstimator* kalmanEstimator, MatrixXd preamble, uint smoothingLag, LinearDetector *linearDetector, std::vector<double> ARcoefficients, bool substractContributionFromKnownSymbols = false);
};

#endif // LINEARFILTERKFBASEDALGORITHM_H
