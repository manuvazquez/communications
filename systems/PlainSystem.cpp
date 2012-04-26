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


#include "PlainSystem.h"

PlainSystem::PlainSystem()
	: BaseSystem()
{
	_powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);

	_kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance);
	_MMSEdetector = new MMSEDetector(_L*(_d+1),_N*(_d+1),_alphabet->variance(),_N*(_d+1));
	_decoratedKalmanEstimator = new LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator(_kalmanEstimator,_MMSEdetector,_alphabet->variance());
}

PlainSystem::~PlainSystem()
{
	delete _kalmanEstimator;
	delete _MMSEdetector;
	delete _decoratedKalmanEstimator;
	
	delete _powerProfile;
}

void PlainSystem::addAlgorithms()
{
	_algorithms.push_back(new LinearFilterBasedAlgorithm("Kalman Filter with noise variance",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients[0]));
	_algorithms.push_back(new LinkedKalmanFilterAndLinearFilterBasedAlgorithm("Kalman Filter with ADJUSTED variance",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_decoratedKalmanEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients[0]));
}

