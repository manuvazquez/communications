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
    _adjustParticlesNumberFromSurvivors = true;
	_adjustSurvivorsFromParticlesNumber = false;

    _powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);
    
    adjustParticlesSurvivors(_nParticles,nSurvivors,_adjustParticlesNumberFromSurvivors,_adjustSurvivorsFromParticlesNumber);

    kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,_ARcoefficients,_ARvariance);

    ResamplingCriterion resamplingCriterion(_resamplingRatio);
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

	_algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (deterministic)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,_nParticles,bestParticlesResamplingAlgorithm,_powerProfile->means(),_powerProfile->variances()/*,ARcoefficients[0]*/));

	_algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (without replacement resampling)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,_nParticles,withoutReplacementResamplingAlgorithm,_powerProfile->means(),_powerProfile->variances()/*,ARcoefficients[0]*/));

	_algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (residual resampling)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,_nParticles,_resamplingAlgorithm,_powerProfile->means(),_powerProfile->variances()/*,ARcoefficients[0]*/));

    _algorithms.push_back(new ViterbiAlgorithm("Viterbi",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (_channel)),_preamble,_d));
}

void PSPvsPSPBasedSMCSystem::saveFrameResults()
{
    SMCSystem::saveFrameResults();
    Octave::toOctaveFileStream(nSurvivors,"nSurvivors",_f);
}

