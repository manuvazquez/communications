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
    nSurvivors = 2;
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
	{
		RLSchannelEstimators.push_back(new RLSEstimator(channelOrderCoefficientsMeans[iChannelOrder],N,forgettingFactor));
		kalmanChannelEstimators.push_back(new KalmanEstimator(channelOrderCoefficientsMeans[iChannelOrder],channelOrderCoefficientsVariances[iChannelOrder],N,ARcoefficients[0],ARvariance));
	}

    ResamplingCriterion resamplingCriterion(resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
	bestParticlesResamplingAlgorithm = new BestParticlesResamplingAlgorithm(resamplingCriterion);

	kalmanEstimatedChannel = estimatedChannel = subestimatedChannel = overestimatedChannel = NULL;

    kalmanEstimator = new KalmanEstimator(powerProfile->Means(),powerProfile->Variances(),N,ARcoefficients[0],ARvariance);
	knownChannelChannelMatrixEstimator = NULL;
	knownChannelChannelMatrixEstimatorEstimatedChannel = NULL;
}


PSPBasedChannelOrderEstimationSystem::~PSPBasedChannelOrderEstimationSystem()
{
// 	delete channel;
	delete powerProfile;

	delete rmmseDetector;

	delete rlsEstimator;
	for(uint iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		delete RLSchannelEstimators[iChannelOrder];
		delete kalmanChannelEstimators[iChannelOrder];
	}

	delete withoutReplacementResamplingAlgorithm;
	delete bestParticlesResamplingAlgorithm;

	delete estimatedChannel;
	delete subestimatedChannel;
	delete overestimatedChannel;

	delete kalmanEstimatedChannel;
	delete kalmanEstimator;

	delete knownChannelChannelMatrixEstimator;
	delete knownChannelChannelMatrixEstimatorEstimatedChannel;
}

void PSPBasedChannelOrderEstimationSystem::BuildChannel()
{
//     channel = new ARchannel(N,L,m,symbols.cols(),ARprocess(powerProfile->GenerateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
// 	channel = new BesselChannel(N,L,m,symbols.cols(),50,2e9,1.0/500.0e3,*powerProfile);
	channel = new TimeInvariantChannel(N,L,m,symbols.cols(),powerProfile->GenerateChannelMatrix(randomGenerator));
}

void PSPBasedChannelOrderEstimationSystem::AddAlgorithms()
{
	ChannelOrderEstimationSystem::AddAlgorithms();

	delete kalmanEstimatedChannel;
	kalmanEstimatedChannel = new EstimatedMIMOChannel(N,L,m,symbols.cols(),preambleLength,kalmanEstimator,symbols,observaciones,ruido->Variances());

	delete knownChannelChannelMatrixEstimator;
	knownChannelChannelMatrixEstimator = new KnownChannelChannelMatrixEstimator(*channel,preambleLength,N);

	delete knownChannelChannelMatrixEstimatorEstimatedChannel;
	knownChannelChannelMatrixEstimatorEstimatedChannel = new EstimatedMIMOChannel(N,L,m,symbols.cols(),preambleLength,knownChannelChannelMatrixEstimator,symbols,observaciones,ruido->Variances());

	delete estimatedChannel;
	estimatedChannel = new EstimatedMIMOChannel(N,L,m,symbols.cols(),preambleLength,rlsEstimator,symbols,observaciones,ruido->Variances());

	delete subestimatedChannel;
	subestimatedChannel = new EstimatedMIMOChannel(N,L,candidateChannelOrders[iTrueChannelOrder-1],symbols.cols(),preambleLength,RLSchannelEstimators[iTrueChannelOrder-1],symbols,observaciones,ruido->Variances());

	delete overestimatedChannel;
	overestimatedChannel = new EstimatedMIMOChannel(N,L,candidateChannelOrders[iTrueChannelOrder+1],symbols.cols(),preambleLength,RLSchannelEstimators[iTrueChannelOrder+1],symbols,observaciones,ruido->Variances());

	algorithms.push_back(new CMEBasedAlgorithm("CME based algorithm",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,preamble,preamble.cols(),symbols));

//     algorithms.push_back(new LinearFilterBasedSMCAlgorithm("Linear Filter Based SMC Algorithm (RLS + MMSE)",*alphabet,L,N,lastSymbolVectorInstant,m,rlsEstimator,rmmseDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

	algorithms.push_back(new UTrellisSearchAlgorithm("UTrellisSearchAlgorithm",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,preamble,preamble.cols(),d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

// 	algorithms.push_back(new UTSAlgorithm("UTSAlgorithm",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,preamble,preamble.cols(),d,nParticles,withoutReplacementResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

// 	algorithms.push_back(new UTSFeedBackAlgorithm("UTSFeedBackAlgorithm (deterministic)",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,preamble,preamble.cols(),d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

// 	algorithms.push_back(new UTSAlgorithm("UTSAlgorithm (deterministic)",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,preamble,preamble.cols(),d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

// 	algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*alphabet,L,N,lastSymbolVectorInstant,m,rlsEstimator,preamble,d,lastSymbolVectorInstant+d,ARcoefficients[0],nSurvivors));

// 	algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,d,nParticles,bestParticlesResamplingAlgorithm,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0]));

    algorithms.push_back(new ViterbiAlgorithm("Viterbi",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d));

// 	algorithms.push_back(new ViterbiAlgorithm("Viterbi (estimated channel)",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (estimatedChannel)),preamble,d));

// 	algorithms.push_back(new ViterbiAlgorithm("Viterbi (Kalman estimated channel)",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (kalmanEstimatedChannel)),preamble,d));

// 	algorithms.push_back(new ViterbiAlgorithm("Viterbi (subestimated channel)",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (subestimatedChannel)),preamble,candidateChannelOrders[iTrueChannelOrder-1]-1));

// 	algorithms.push_back(new ViterbiAlgorithm("Viterbi (overestimated channel)",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (overestimatedChannel)),preamble,candidateChannelOrders[iTrueChannelOrder+1]-1));

// 	algorithms.push_back(new ISIS("ISIS",*alphabet,L,N,lastSymbolVectorInstant,kalmanChannelEstimators,preamble,preamble.cols(),d,nParticles,algoritmoRemuestreo));
}

void PSPBasedChannelOrderEstimationSystem::BeforeEndingFrame(int iFrame)
{
    ChannelOrderEstimationSystem::BeforeEndingFrame(iFrame);
    Util::ScalarToOctaveFileStream(nSurvivors,"nSurvivors",f);
	Util::ScalarToOctaveFileStream(forgettingFactor,"forgettingFactor",f);
	Util::ScalarToOctaveFileStream(forgettingFactorDetector,"forgettingFactorDetector",f);
}

