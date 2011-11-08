/***************************************************************************
 *   Copyright (C) 2006 by Manu   *
 *   manu@rustneversleeps   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "PSPvsPSPBasedSMCSystem.h"

PSPvsPSPBasedSMCSystem::PSPvsPSPBasedSMCSystem()
 : SMCSystem()
{
    nSurvivors = 2;
    adjustParticlesNumberFromSurvivors = true;

    _powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);

	if(adjustParticlesNumberFromSurvivors)
	{
		nParticles = (int)pow((double)_alphabet->length(),_N*(_m-1))*nSurvivors;
        cout << "Number of particles adjusted to " << nParticles << endl;
    }

    kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,ARcoefficients,ARvariance);

    ResamplingCriterion resamplingCriterion(resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
	bestParticlesResamplingAlgorithm = new BestParticlesResamplingAlgorithm(resamplingCriterion);
}


PSPvsPSPBasedSMCSystem::~PSPvsPSPBasedSMCSystem()
{
	delete _powerProfile;
	delete kalmanEstimator;

	delete withoutReplacementResamplingAlgorithm;
	delete bestParticlesResamplingAlgorithm;
}

void PSPvsPSPBasedSMCSystem::addAlgorithms()
{
	_algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,_iLastSymbolVectorToBeDetected+_d,nSurvivors));

	_algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (deterministic)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,nParticles,bestParticlesResamplingAlgorithm,_powerProfile->means(),_powerProfile->variances()/*,ARcoefficients[0]*/));

	_algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (without replacement resampling)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,nParticles,withoutReplacementResamplingAlgorithm,_powerProfile->means(),_powerProfile->variances()/*,ARcoefficients[0]*/));

	_algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (residual resampling)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,nParticles,algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances()/*,ARcoefficients[0]*/));

    _algorithms.push_back(new ViterbiAlgorithm("Viterbi",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (_channel)),_preamble,_d));
}

void PSPvsPSPBasedSMCSystem::buildSystemSpecificVariables()
{
//     channel = new ARchannel(N,L,m,symbols.cols(),ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
	_channel = new BesselChannel(_N,_L,_m,_symbols.cols(),50,2e9,1.0/500.0e3,*_powerProfile);
}

void PSPvsPSPBasedSMCSystem::beforeEndingFrame()
{
    SMCSystem::beforeEndingFrame();
    Util::scalarToOctaveFileStream(nSurvivors,"nSurvivors",_f);
}

