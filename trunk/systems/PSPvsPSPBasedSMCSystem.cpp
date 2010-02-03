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

    powerProfile = new FlatPowerProfile(L,N,m,1.0);

	if(adjustParticlesNumberFromSurvivors)
	{
		nParticles = (int)pow((double)alphabet->length(),N*(m-1))*nSurvivors;
        cout << "Number of particles adjusted to " << nParticles << endl;
    }

    kalmanEstimator = new KalmanEstimator(powerProfile->means(),powerProfile->variances(),N,ARcoefficients,ARvariance);

    ResamplingCriterion resamplingCriterion(resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
	bestParticlesResamplingAlgorithm = new BestParticlesResamplingAlgorithm(resamplingCriterion);
}


PSPvsPSPBasedSMCSystem::~PSPvsPSPBasedSMCSystem()
{
	delete powerProfile;
	delete kalmanEstimator;

	delete withoutReplacementResamplingAlgorithm;
	delete bestParticlesResamplingAlgorithm;
}

void PSPvsPSPBasedSMCSystem::AddAlgorithms()
{
	algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,iLastSymbolVectorToBeDetected+d,ARcoefficients[0],nSurvivors));

	algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (deterministic)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,nParticles,bestParticlesResamplingAlgorithm,powerProfile->means(),powerProfile->variances(),ARcoefficients[0]));

	algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (without replacement resampling)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,nParticles,withoutReplacementResamplingAlgorithm,powerProfile->means(),powerProfile->variances(),ARcoefficients[0]));

	algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (residual resampling)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),ARcoefficients[0]));

    algorithms.push_back(new ViterbiAlgorithm("Viterbi",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d));
}

void PSPvsPSPBasedSMCSystem::BuildChannel()
{
//     channel = new ARchannel(N,L,m,symbols.cols(),ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
	channel = new BesselChannel(N,L,m,symbols.cols(),50,2e9,1.0/500.0e3,*powerProfile);
}

void PSPvsPSPBasedSMCSystem::BeforeEndingFrame(int iFrame)
{
    SMCSystem::BeforeEndingFrame(iFrame);
    Util::scalarToOctaveFileStream(nSurvivors,"nSurvivors",f);
}

