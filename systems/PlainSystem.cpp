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
// 	_ARcoefficients = ARprocess::parametersFromYuleWalker(2,_velocity,_carrierFrequency,_T,_ARvariance);
// 	
// 	std::cout << "AR coefficients computed from Yule-Walker equations:" << std::endl << _ARcoefficients << std::endl;
// 	std::cout << "Variance = " << _ARvariance << std::endl;
	
	_powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);
	
// 	_MMSEdetector = new MMSEDetector(_L*(_d+1),_N*(_d+1),_alphabet->variance(),_N*(_d+1));
	_MMSEdetector = new MMSEDetector(_L*(_d+1),_N*(_m+_d),_alphabet->variance(),_N*(_d+1));

	_kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance);
	_linearFilterAwareNoiseVarianceAdjustingKalmanEstimator = new LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance,_MMSEdetector,_alphabet->variance());
	
// 	_kalmanFilterAwareMMSEDetector = new KalmanFilterAwareMMSEDetector(_L*(_d+1),_N*(_d+1),_alphabet->variance(),_N*(_d+1),_kalmanEstimator);
	_kalmanFilterAwareMMSEDetector = new KalmanFilterAwareMMSEDetector(_L*(_d+1),_N*(_m+_d),_alphabet->variance(),_N*(_d+1),_kalmanEstimator,_ARcoefficients);
	
	_knownChannelChannelMatrixEstimator = NULL;
}

PlainSystem::~PlainSystem()
{
	delete _kalmanEstimator;
	delete _linearFilterAwareNoiseVarianceAdjustingKalmanEstimator;
	
	delete _MMSEdetector;
	delete _knownChannelChannelMatrixEstimator;
	
	delete _kalmanFilterAwareMMSEDetector;
	
	delete _powerProfile;
}

void PlainSystem::addAlgorithms()
{
	_algorithms.push_back(new KalmanFilterAwareMMSEBasedAlgorithm("Kalman Filter-aware-MMSE + Kalman Filter",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_kalmanFilterAwareMMSEDetector,_ARcoefficients));
	_algorithms.push_back(new LinearFilterBasedAlgorithm("MMSE + Kalman Filter",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));

// 	_algorithms.push_back(new LinkedKalmanFilterAndLinearFilterBasedAlgorithm("Linear Filter + Kalman Filter with ADJUSTED variance",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_linearFilterAwareNoiseVarianceAdjustingKalmanEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));

	_algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_symbols));
	
    delete _knownChannelChannelMatrixEstimator;
    _knownChannelChannelMatrixEstimator = new KnownChannelChannelMatrixEstimator(_channel,_preambleLength,_N);
	_algorithms.push_back(new LinearFilterBasedAlgorithm("MMSE with Known Channel",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownChannelChannelMatrixEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));
}
