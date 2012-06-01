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
	
	double yuleWalkerARvariance;
	std::vector<double> yuleWalkerARcoefficients = ARprocess::parametersFromYuleWalker(_ARcoefficients.size(),_velocity,_carrierFrequency,_T,yuleWalkerARvariance);

	std::cout << "AR coefficients computed from Yule-Walker equations:" << std::endl << yuleWalkerARcoefficients << std::endl;
	std::cout << "Variance = " << yuleWalkerARvariance << std::endl;
	
// 	_ARcoefficients = yuleWalkerARcoefficients;
// 	_ARvariance = yuleWalkerARvariance;
	
	_powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);
	
	_MMSEdetector = new MMSEDetector(_L*(_d+1),_N*(_m+_d),_alphabet->variance(),_N*(_d+1));

	_kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance);
	
	_SERawareKalmanEstimatorDecorator = new SERawareKalmanEstimatorDecorator(_kalmanEstimator,noiseVariances,SERs,_alphabet->differencesBetweenSymbols());
	
	_kalmanFilterAwareMMSEDetector = new KalmanFilterAwareMMSEDetector(_L*(_d+1),_N*(_m+_d),_alphabet->variance(),_N*(_d+1),_kalmanEstimator,_ARcoefficients);
	
	_knownChannelChannelMatrixEstimator = NULL;
	
	_knownSymbolsKalmanEstimator = NULL;
}

PlainSystem::~PlainSystem()
{
	delete _kalmanEstimator;
	delete _knownSymbolsKalmanEstimator;
	delete _SERawareKalmanEstimatorDecorator;
	
	delete _MMSEdetector;
	delete _knownChannelChannelMatrixEstimator;
	
	delete _kalmanFilterAwareMMSEDetector;
	
	delete _powerProfile;
}

void PlainSystem::addAlgorithms()
{
	_algorithms.push_back(new LinearFilterBasedAlgorithm("MMSE + Kalman Filter",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));
	_algorithms.push_back(new LinearFilterBasedAlgorithm("MMSE + Known SER Kalman Filter",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_SERawareKalmanEstimatorDecorator,_preamble,_d,_MMSEdetector,_ARcoefficients));
	_algorithms.push_back(new KalmanFilterAwareMMSEBasedAlgorithm("KF-aware MMSE + Kalman Filter",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_kalmanFilterAwareMMSEDetector,_ARcoefficients));
	_algorithms.push_back(new KalmanFilterAwareMMSEBasedAlgorithm("KF-aware MMSE + Known SER Kalman Filter",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_SERawareKalmanEstimatorDecorator,_preamble,_d,_kalmanFilterAwareMMSEDetector,_ARcoefficients));
	
	delete _knownSymbolsKalmanEstimator;
	_knownSymbolsKalmanEstimator = new KnownSymbolsKalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance,_symbols,_preambleLength);
	_algorithms.push_back(new LinearFilterBasedAlgorithm("MMSE + Known Symbols Kalman Filter",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownSymbolsKalmanEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));
	_algorithms.push_back(new KalmanFilterAwareMMSEBasedAlgorithm("KF-aware MMSE + Known Symbols Kalman Filter",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownSymbolsKalmanEstimator,_preamble,_d,_kalmanFilterAwareMMSEDetector,_ARcoefficients));
	
    delete _knownChannelChannelMatrixEstimator;
    _knownChannelChannelMatrixEstimator = new KnownChannelChannelMatrixEstimator(_channel,_preambleLength,_N);
	_algorithms.push_back(new LinearFilterBasedAlgorithm("MMSE with Known Channel at Detection Time",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownChannelChannelMatrixEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));
	_algorithms.push_back(new LinearFilterBasedAlgorithmWithKnownChannel("MMSE with Known Channel at Every Time Instant",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_knownChannelChannelMatrixEstimator,_preamble,_d,_MMSEdetector,_ARcoefficients));
	
	_algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_symbols));
}
