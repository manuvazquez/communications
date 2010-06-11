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

#include <bashcolors.h>

WSA08System::WSA08System()
 : ChannelOrderEstimationSystem()
{
    nSurvivors = 1;
    adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;

    forgettingFactor = 0.99;
    forgettingFactorDetector = 0.95;

	velocity = 50.0; // m/s

    _powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);
// 	powerProfile = new ExponentialPowerProfile(L,N,m,1.8e-6,1.0/500.0e3);

	if(adjustParticlesNumberFromSurvivors)
	{
		nParticles = (int)pow((double)_alphabet->length(),_N*(_m-1))*nSurvivors;
        cout << COLOR_WHITE << "Number of particles adjusted to " << COLOR_NORMAL << nParticles << endl;
    }

	if(adjustSurvivorsFromParticlesNumber)
	{
		cout << COLOR_WHITE << "Number of survivors adjusted from " << COLOR_NORMAL << nSurvivors;
		nSurvivors = int(ceil(double(nParticles)/pow(2.0,double(_N*(_m-1)))));
		cout << COLOR_WHITE << " to " << COLOR_NORMAL << nSurvivors << endl;
	}

	rmmseDetector = new RMMSEDetector(_L*(c+_d+1),_N*(_m+c+_d),_alphabet->variance(),forgettingFactorDetector,_N*(_d+1));

	rlsEstimator = new RLSEstimator(_powerProfile->means(),_N,forgettingFactor);
	for(uint iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		RLSchannelEstimators.push_back(new RLSEstimator(channelOrderCoefficientsMeans[iChannelOrder],_N,forgettingFactor));
		kalmanChannelEstimators.push_back(new KalmanEstimator(channelOrderCoefficientsMeans[iChannelOrder],channelOrderCoefficientsVariances[iChannelOrder],_N,ARcoefficients,ARvariance));
		noForgetRLSchannelEstimators.push_back(new RLSEstimator(channelOrderCoefficientsMeans[iChannelOrder],_N,1.0));
	}

    ResamplingCriterion resamplingCriterion(resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
	bestParticlesResamplingAlgorithm = new BestParticlesResamplingAlgorithm(resamplingCriterion);

    kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,ARcoefficients,ARvariance);
}


WSA08System::~WSA08System()
{
// 	delete channel;
	delete _powerProfile;

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

void WSA08System::buildChannel()
{
//     channel = new ARchannel(N,L,m,symbols.cols(),ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
	_channel = new BesselChannel(_N,_L,_m,_symbols.cols(),velocity,2e9,1.0/500.0e3,*_powerProfile);
	
	
// 	std::vector<uint> subchannelOrders(3,1);
// 	subchannelOrders[0] = 3;
// 	subchannelOrders[2] = 1;
// 	
// 	dynamic_cast<StillMemoryMIMOChannel*>(_channel)->setSubchannelOrders(subchannelOrders);
// 	cout << _channel->at(100) << endl;
// 	getchar();

// 	channel = new TimeInvariantChannel(N,L,m,symbols.cols(),powerProfile->generateChannelMatrix(randomGenerator));
#ifdef DEBUG
	cout << "El canal al principio" << endl << (*channel)[preamble.cols()];
	cout << "El canal al final" << endl << (*channel)[frameLength];
#endif
}

void WSA08System::addAlgorithms()
{
	ChannelOrderEstimationSystem::addAlgorithms();

// 	algorithms.push_back(new MLSDmAlgorithm("MLSDmAlgorithm",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,RLSchannelEstimators,preamble,preamble.cols(),d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

	_algorithms.push_back(new MLSDmAlgorithm("MLSD-m",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,kalmanChannelEstimators,_preamble,_preamble.cols(),_d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

	_algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,_iLastSymbolVectorToBeDetected+_d,nSurvivors));
	
// 	algorithms.push_back(new PSPAlgorithm("PSPAlgorithm +1 survivor",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,iLastSymbolVectorToBeDetected+d,nSurvivors+1));

// 	algorithms.push_back(new PSPAlgorithm("PSPAlgorithm (underestimated)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,candidateChannelOrders[iTrueChannelOrder-1],dynamic_cast <KalmanEstimator *> (kalmanChannelEstimators[iTrueChannelOrder-1]),preamble,candidateChannelOrders[iTrueChannelOrder-1]-1,iLastSymbolVectorToBeDetected+candidateChannelOrders[iTrueChannelOrder-1]-1,nSurvivors));

// 	algorithms.push_back(new PSPAlgorithm("PSPAlgorithm (overestimated)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,candidateChannelOrders[iTrueChannelOrder+1],dynamic_cast <KalmanEstimator *> (kalmanChannelEstimators[iTrueChannelOrder+1]),preamble,candidateChannelOrders[iTrueChannelOrder+1]-1,iLastSymbolVectorToBeDetected+candidateChannelOrders[iTrueChannelOrder+1]-1,nSurvivors));

    _algorithms.push_back(new ViterbiAlgorithm("Viterbi",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (_channel)),_preamble,_d));
}

void WSA08System::beforeEndingFrame()
{
    ChannelOrderEstimationSystem::beforeEndingFrame();
    Util::scalarToOctaveFileStream(nSurvivors,"nSurvivors",_f);
	Util::scalarToOctaveFileStream(forgettingFactor,"forgettingFactor",_f);
	Util::scalarToOctaveFileStream(forgettingFactorDetector,"forgettingFactorDetector",_f);
	Util::scalarToOctaveFileStream(velocity,"velocity",_f);
}

