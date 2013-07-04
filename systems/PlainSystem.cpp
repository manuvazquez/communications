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
	std::vector<double> noiseVariances;
	std::vector<double> SERs;
	
	// -------------------------------- parameters
	
	xml_node<> *thisSystemParameters = get_child(_doc.first_node(),"PlainSystem");
	
	if(!thisSystemParameters)
		throw RuntimeException("PlainSystem::PlainSystem: cannot find parameters for this system.");
	
	readMultiValuedParameterFromXML(thisSystemParameters,"noiseVariances",noiseVariances);
	readMultiValuedParameterFromXML(thisSystemParameters,"SERs",SERs);
	
	// --------------------------------
	
	_powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);
	
	_MMSEdetector = new MMSEDetector(_L*(_d+1),_N*(_m+_d),_alphabet->variance(),_N*(_d+1));
	
    _ICMMSEdetector = new MMSEDetector(_L*(_d+1),_N*(_d+1),_alphabet->variance(),_N*(_d+1));

	_kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance);
	
	_kalmanFilterAwareMMSEDetector = new KalmanFilterAwareMMSEDetector(_L*(_d+1),_N*(_m+_d),_alphabet->variance(),_N*(_d+1),_kalmanEstimator,_ARcoefficients);
	
	_ICKFAwareMMSEDetector = new KalmanFilterAwareMMSEDetector(_L*(_d+1),_N*(_d+1),_alphabet->variance(),_N*(_d+1),_kalmanEstimator,_ARcoefficients,true);
	
	_knownChannelChannelMatrixEstimator = NULL;
	
	_knownSymbolsKalmanEstimator = NULL;
}

PlainSystem::~PlainSystem()
{
	delete _kalmanEstimator;
	delete _knownSymbolsKalmanEstimator;
	
	delete _MMSEdetector;
	delete _knownChannelChannelMatrixEstimator;
	
	delete _kalmanFilterAwareMMSEDetector;
	
	delete _ICKFAwareMMSEDetector;
	
	delete _powerProfile;
}

void PlainSystem::addAlgorithms()
{
	_algorithms.push_back(new LinearFilterKFBasedAlgorithm("MMSE + KF",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));
	_algorithms.push_back(new LinearFilterKFBasedAlgorithm("MMSE + KF with IC",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_ICMMSEdetector,_ARcoefficients,true));
	_algorithms.push_back(new KalmanFilterAwareMMSEBasedAlgorithm("SOS-MMSE + KF with IC",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_ICKFAwareMMSEDetector,_ARcoefficients,true));
	_algorithms.push_back(new KalmanFilterAwareMMSEBasedAlgorithm("SOS-MMSE + KF",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_kalmanFilterAwareMMSEDetector,_ARcoefficients));
	
	delete _knownSymbolsKalmanEstimator;
	_knownSymbolsKalmanEstimator = new KnownSymbolsKalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance,_symbols,_preambleLength);
	
	_algorithms.push_back(new LinearFilterKFBasedAlgorithm("MMSE + Known Symbols KF",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownSymbolsKalmanEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));
	_algorithms.push_back(new KalmanFilterAwareMMSEBasedAlgorithm("SOS-MMSE + Known Symbols KF",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownSymbolsKalmanEstimator,_preamble,_d,_kalmanFilterAwareMMSEDetector,_ARcoefficients));
	
    delete _knownChannelChannelMatrixEstimator;
    _knownChannelChannelMatrixEstimator = new KnownChannelChannelMatrixEstimator(_channel,_preambleLength,_N);
	
// 	_algorithms.push_back(new LinearFilterBasedAlgorithm("MMSE with Known Channel at Detection Time",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownChannelChannelMatrixEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));
	_algorithms.push_back(new LinearFilterBasedAlgorithm("MMSE with Known Channel at Detection Time",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownChannelChannelMatrixEstimator,_preamble,_d,_ICMMSEdetector,_ARcoefficients,true));
	
	_algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_symbols));
}
