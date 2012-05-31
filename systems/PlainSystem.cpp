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
	double yuleWalkerARvariance;
	std::vector<double> yuleWalkerARcoefficients = ARprocess::parametersFromYuleWalker(_ARcoefficients.size(),_velocity,_carrierFrequency,_T,yuleWalkerARvariance);

	std::cout << "AR coefficients computed from Yule-Walker equations:" << std::endl << yuleWalkerARcoefficients << std::endl;
	std::cout << "Variance = " << yuleWalkerARvariance << std::endl;
	
// 	_ARcoefficients = yuleWalkerARcoefficients;
// 	_ARvariance = yuleWalkerARvariance;
	
	_powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);
	
	_MMSEdetector = new MMSEDetector(_L*(_d+1),_N*(_m+_d),_alphabet->variance(),_N*(_d+1));

	_kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance);
	
	// for BPSK...and from previous simulation results
	std::vector<double> noiseVariances;
	noiseVariances.push_back(3.00712);noiseVariances.push_back(1.50713);noiseVariances.push_back(0.755355);
	noiseVariances.push_back(0.378574);noiseVariances.push_back(0.189737);noiseVariances.push_back(0.0950936);
	noiseVariances.push_back(0.0476597);noiseVariances.push_back(0.0238864);noiseVariances.push_back(0.0119716);
	noiseVariances.push_back(0.006);noiseVariances.push_back(0.00300712);noiseVariances.push_back(0.00150713);
	noiseVariances.push_back(0.000755355);
	std::vector<double> SERs;
	SERs.push_back(0.0982456);SERs.push_back(0.0596491);SERs.push_back(0.0362573);SERs.push_back(0.0245614);
	SERs.push_back(0.0175439);SERs.push_back(0.0140351);SERs.push_back(0.0116959);SERs.push_back(0.00584795);
	SERs.push_back(0.00350877);SERs.push_back(0.00233918);SERs.push_back(0.00233918);SERs.push_back(0.00233918);
	SERs.push_back(0.00233918);
	
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
	_algorithms.push_back(new KalmanFilterAwareMMSEBasedAlgorithm("KF-aware MMSE + Kalman Filter",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimator,_preamble,_d,_kalmanFilterAwareMMSEDetector,_ARcoefficients));
	
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
