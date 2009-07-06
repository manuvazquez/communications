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
#include "WSA08System.h"

WSA08System::WSA08System()
 : ChannelOrderEstimationSystem()
{
    nSurvivors = 3;
	adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;

    forgettingFactor = 0.99;
    forgettingFactorDetector = 0.95;

	velocity = 50.0; // m/s


    powerProfile = new FlatPowerProfile(L,N,m,1.0);
// 	powerProfile = new ExponentialPowerProfile(L,N,m,1.8e-6,1.0/500.0e3);

	if(adjustParticlesNumberFromSurvivors)
	{
		nParticles = (int)pow((double)alphabet->length(),N*(m-1))*nSurvivors;
        cout << "Number of particles adjusted to " << nParticles << endl;
    }

	if(adjustSurvivorsFromParticlesNumber)
	{
		cout << "number of survivors adjusted from " << nSurvivors;
		nSurvivors = int(ceil(double(nParticles)/pow(2.0,double(N*(m-1)))));
		cout << " to " << nSurvivors << endl;
	}

	rmmseDetector = new RMMSEDetector(L*(c+d+1),N*(m+c+d),alphabet->variance(),forgettingFactorDetector,N*(d+1));

	rlsEstimator = new RLSEstimator(powerProfile->means(),N,forgettingFactor);
	for(uint iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		RLSchannelEstimators.push_back(new RLSEstimator(channelOrderCoefficientsMeans[iChannelOrder],N,forgettingFactor));
		kalmanChannelEstimators.push_back(new KalmanEstimator(channelOrderCoefficientsMeans[iChannelOrder],channelOrderCoefficientsVariances[iChannelOrder],N,ARcoefficients,ARvariance));
		noForgetRLSchannelEstimators.push_back(new RLSEstimator(channelOrderCoefficientsMeans[iChannelOrder],N,1.0));
	}

    ResamplingCriterion resamplingCriterion(resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
	bestParticlesResamplingAlgorithm = new BestParticlesResamplingAlgorithm(resamplingCriterion);

    kalmanEstimator = new KalmanEstimator(powerProfile->means(),powerProfile->variances(),N,ARcoefficients,ARvariance);
}


WSA08System::~WSA08System()
{
// 	delete channel;
	delete powerProfile;

	delete rmmseDetector;

	delete rlsEstimator;
	for(uint iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		delete RLSchannelEstimators[iChannelOrder];
		delete kalmanChannelEstimators[iChannelOrder];
		delete noForgetRLSchannelEstimators[iChannelOrder];
	}

	delete withoutReplacementResamplingAlgorithm;
	delete bestParticlesResamplingAlgorithm;

	delete kalmanEstimator;
}

void WSA08System::BuildChannel()
{
//     channel = new ARchannel(N,L,m,symbols.cols(),ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
	channel = new BesselChannel(N,L,m,symbols.cols(),velocity,2e9,1.0/500.0e3,*powerProfile);
// 	channel = new TimeInvariantChannel(N,L,m,symbols.cols(),powerProfile->generateChannelMatrix(randomGenerator));
#ifdef DEBUG
	cout << "El canal al principio" << endl << (*channel)[preamble.cols()];
	cout << "El canal al final" << endl << (*channel)[frameLength];
#endif
}

void WSA08System::AddAlgorithms()
{
	ChannelOrderEstimationSystem::AddAlgorithms();

// 	algorithms.push_back(new MLSDmAlgorithm("MLSDmAlgorithm",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,RLSchannelEstimators,preamble,preamble.cols(),d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

	algorithms.push_back(new MLSDmAlgorithm("MLSD-m",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,kalmanChannelEstimators,preamble,preamble.cols(),d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

	algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,iLastSymbolVectorToBeDetected+d,ARcoefficients[0],nSurvivors));

	algorithms.push_back(new PSPAlgorithm("PSPAlgorithm (underestimated)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,candidateChannelOrders[iTrueChannelOrder-1],kalmanChannelEstimators[iTrueChannelOrder-1],preamble,candidateChannelOrders[iTrueChannelOrder-1]-1,iLastSymbolVectorToBeDetected+candidateChannelOrders[iTrueChannelOrder-1]-1,ARcoefficients[0],nSurvivors));

	algorithms.push_back(new PSPAlgorithm("PSPAlgorithm (overestimated)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,candidateChannelOrders[iTrueChannelOrder+1],kalmanChannelEstimators[iTrueChannelOrder+1],preamble,candidateChannelOrders[iTrueChannelOrder+1]-1,iLastSymbolVectorToBeDetected+candidateChannelOrders[iTrueChannelOrder+1]-1,ARcoefficients[0],nSurvivors));

    algorithms.push_back(new ViterbiAlgorithm("Viterbi",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d));
}

void WSA08System::BeforeEndingFrame(int iFrame)
{
    ChannelOrderEstimationSystem::BeforeEndingFrame(iFrame);
    Util::ScalarToOctaveFileStream(nSurvivors,"nSurvivors",f);
	Util::ScalarToOctaveFileStream(forgettingFactor,"forgettingFactor",f);
	Util::ScalarToOctaveFileStream(forgettingFactorDetector,"forgettingFactorDetector",f);
	Util::ScalarToOctaveFileStream(velocity,"velocity",f);
}

