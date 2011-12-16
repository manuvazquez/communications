/*
    Copyright 2011 Manu <manuavazquez@gmail.com>

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


#include "ARparametersComparisonCDMASystem.h"

#include <defines.h>
#include <OldUnknownActiveUsersLinearFilterBasedSMCAlgorithm.h>

void ARparametersComparisonCDMASystem::addAlgorithms()
{	
	for(uint i=0;i<_ARCoeffs.size();i++)
	{	
		delete _cdmaKalmanEstimator; 
		_cdmaKalmanEstimator = new CDMAKalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_ARCoeffs[i],_ARvariances[i],_spreadingCodes);

// 		std::stringstream algorithmName;	
// 		algorithmName << "KF coeff[0]=" << _ARCoeffs[i][0] << " coeff[1]=" << _ARCoeffs[i][1] << " var=" << _ARvariances[i];
// 		_algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm(algorithmName.str(),*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_preamble,_symbols));
		
		std::stringstream algorithmName;	
		algorithmName << "SMC-LF coeff[0]=" << _ARCoeffs[i][0] << " coeff[1]=" << _ARCoeffs[i][1] << " var=" << _ARvariances[i];
		_algorithms.push_back(new OldUnknownActiveUsersLinearFilterBasedSMCAlgorithm (algorithmName.str(),*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_mmseDetector,_preamble,_d,nParticles,algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_usersActivityPdfs));
	}
}

ARparametersComparisonCDMASystem::ARparametersComparisonCDMASystem()
{
	std::vector<double> ARcoeffs(2);
	
	ARcoeffs[0] = 0.59999;
	ARcoeffs[1] = 0.39999;
	_ARvariances.push_back(0.0001);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.59999;
	ARcoeffs[1] = 0.39999;
	_ARvariances.push_back(0.00001);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.59999;
	ARcoeffs[1] = 0.39999;
	_ARvariances.push_back(0.000001);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.59999;
	ARcoeffs[1] = 0.39999;
	_ARvariances.push_back(0.0000001);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.59999;
	ARcoeffs[1] = 0.39999;
	_ARvariances.push_back(0.001);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.59999;
	ARcoeffs[1] = 0.39999;
	_ARvariances.push_back(0.01);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.59999;
	ARcoeffs[1] = 0.39999;
	_ARvariances.push_back(0.1);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.69999;
	ARcoeffs[1] = 0.29999;
	_ARvariances.push_back(0.0001);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.49999;
	ARcoeffs[1] = 0.49999;
	_ARvariances.push_back(0.0001);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.39999;
	ARcoeffs[1] = 0.59999;
	_ARvariances.push_back(0.0001);
	_ARCoeffs.push_back(ARcoeffs);
	
	ARcoeffs[0] = 0.39999;
	ARcoeffs[1] = 0.59999;
	_ARvariances.push_back(0.01);
	_ARCoeffs.push_back(ARcoeffs);
	
	_ARvariances.push_back(FUNNY_VALUE);
	ARcoeffs = ARprocess::parametersFromYuleWalker(2,_velocity,_carrierFrequency,_T,_ARvariances[_ARvariances.size()-1]);
	_ARCoeffs.push_back(ARcoeffs);
}
