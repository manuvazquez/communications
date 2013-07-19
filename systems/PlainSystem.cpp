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
	
    _ICMMSEdetector = new MMSEDetector(_L*(_d+1),_N*(_d+1),_alphabet->variance(),_N*(_d+1));

	_kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance);
	
	_SOSMMSEDetector = new SOSMMSEDetector(_L*(_d+1),_N*(_m+_d),_alphabet->variance(),_N*(_d+1),_kalmanEstimator,_ARcoefficients);
	
	_ICSOSMMSEDetector = new SOSMMSEDetector(_L*(_d+1),_N*(_d+1),_alphabet->variance(),_N*(_d+1),_kalmanEstimator,_ARcoefficients,true);
	
	_knownChannelChannelMatrixEstimator = NULL;
	
	_knownSymbolsKalmanEstimator = NULL;
}

PlainSystem::~PlainSystem()
{
	delete _kalmanEstimator;
	delete _knownSymbolsKalmanEstimator;
	
	delete _knownChannelChannelMatrixEstimator;
	
	delete _SOSMMSEDetector;
	
	delete _ICSOSMMSEDetector;
	
	delete _powerProfile;
}

void PlainSystem::addAlgorithms()
{
// 	_algorithms.push_back(new LinearFilterKFBasedAlgorithm("MMSE + KF",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));
	_algorithms.push_back(new LinearFilterKFBasedAlgorithm("MMSE + KF with IC",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_ICMMSEdetector,_ARcoefficients,true));
	_algorithms.push_back(new SOSMMSEBasedAlgorithm("SOS-MMSE + KF with IC",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_ICSOSMMSEDetector,_ARcoefficients,true));
// 	_algorithms.push_back(new KalmanFilterAwareMMSEBasedAlgorithm("SOS-MMSE + KF",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_SOSMMSEDetector,_ARcoefficients));
	
	delete _knownSymbolsKalmanEstimator;
	_knownSymbolsKalmanEstimator = new KnownSymbolsKalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance,_symbols,_preambleLength);
	
	_algorithms.push_back(new LinearFilterKFBasedAlgorithm("MMSE + Known Symbols KF with IC",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownSymbolsKalmanEstimator,_preamble,_d,_ICMMSEdetector,_ARcoefficients,true));
	_algorithms.push_back(new SOSMMSEBasedAlgorithm("SOS-MMSE + Known Symbols KF with IC",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownSymbolsKalmanEstimator,_preamble,_d,_ICSOSMMSEDetector,_ARcoefficients,true));
	
    delete _knownChannelChannelMatrixEstimator;
    _knownChannelChannelMatrixEstimator = new KnownChannelChannelMatrixEstimator(_channel,_preambleLength,_N);
	
	// the channel is only known at detection time (channel matrices corresponding to the subsequent time instants, needed to perform smoothing, are obtained from channel model)
	_algorithms.push_back(new LinearFilterBasedAlgorithm("MMSE + Known Channel with IC",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownChannelChannelMatrixEstimator,_preamble,_d,_ICMMSEdetector,_ARcoefficients,true));
	
	_algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_symbols));
}
