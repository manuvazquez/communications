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
#include "PSPBasedChannelOrderEstimationSystem.h"

PSPBasedChannelOrderEstimationSystem::PSPBasedChannelOrderEstimationSystem()
 : ChannelOrderEstimationSystem()
{
//     nSurvivors = 2;
//     adjustParticlesNumberFromSurvivors = false;

    forgettingFactor = 0.99;
    forgettingFactorDetector = 0.95;


    powerProfile = new FlatPowerProfile(L,N,m,1.0);

// 	if(adjustParticlesNumberFromSurvivors)
// 	{
// 		nParticles = (int)pow((double)alphabet->Length(),N*(m-1))*nSurvivors;
//         cout << "Number of particles adjusted to " << nParticles << endl;
//     }

	rmmseDetector = new RMMSEDetector(L*(c+d+1),N*(m+c+d),alphabet->Variance(),forgettingFactorDetector,N*(d+1));

	rlsEstimator = new RLSEstimator(powerProfile->Means(),N,forgettingFactor);
	for(uint iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
		RLSchannelEstimators.push_back(new RLSEstimator(channelOrderCoefficientsMeans[iChannelOrder],N,forgettingFactor));

    ResamplingCriterion resamplingCriterion(resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
}


PSPBasedChannelOrderEstimationSystem::~PSPBasedChannelOrderEstimationSystem()
{
}

void PSPBasedChannelOrderEstimationSystem::BuildChannel()
{
    channel = new ARchannel(N,L,m,symbols.cols(),ARprocess(powerProfile->GenerateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
}

void PSPBasedChannelOrderEstimationSystem::AddAlgorithms()
{
    algorithms.push_back(new LinearFilterBasedSMCAlgorithm("Linear Filter Based SMC Algorithm (RLS + MMSE)",*alphabet,L,N,lastSymbolVectorInstant,m,rlsEstimator,rmmseDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

	algorithms.push_back(new UPSPBasedSMCAlgorithm("Unknown channel order PSP based SMC algorithm",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,preamble,preamble.cols(),d,nParticles,withoutReplacementResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

// 	algorithms.push_back(new ISIS("ISIS",*alphabet,L,N,lastSymbolVectorInstant,kalmanChannelEstimators,preamble,preamble.cols(),d,nParticles,&algoritmoRemuestreo,canal,symbols));
}

void PSPBasedChannelOrderEstimationSystem::BeforeEndingFrame(int iFrame)
{
    SMCSystem::BeforeEndingFrame(iFrame);
//     Util::ScalarToOctaveFileStream(nSurvivors,"nSurvivors",f);
	Util::ScalarToOctaveFileStream(forgettingFactor,"forgettingFactor",f);
}

