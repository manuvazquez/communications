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
#include "ChannelOrderEstimationSystem.h"

// #define DEBUG

ChannelOrderEstimationSystem::ChannelOrderEstimationSystem()
 : SMCSystem()
{
	candidateChannelOrders.push_back(2);candidateChannelOrders.push_back(3);candidateChannelOrders.push_back(4);
	candidateChannelOrders.push_back(5);
// 	candidateChannelOrders.push_back(6);candidateChannelOrders.push_back(7);

	channelOrderCoefficientsMeans.resize(candidateChannelOrders.size());
	channelOrderCoefficientsVariances.resize(candidateChannelOrders.size());

	iTrueChannelOrder = -1;
	for(uint iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		channelOrderCoefficientsMeans[iChannelOrder] = MatrixXd::Zero(_L,_N*candidateChannelOrders[iChannelOrder]);
		channelOrderCoefficientsVariances[iChannelOrder] = MatrixXd::Ones(_L,_N*candidateChannelOrders[iChannelOrder]);
		if(candidateChannelOrders[iChannelOrder] == _m)
			iTrueChannelOrder = iChannelOrder;
	}

	if(iTrueChannelOrder==-1)
		throw RuntimeException("ChannelOrderEstimationSystem::ChannelOrderEstimationSystem: the memory of the channel is not one of the possible candidates.");

	// channel order APP evolution
    channelOrderAPPsAlongTime.reserve(_nFrames);
}

void ChannelOrderEstimationSystem::beforeEndingFrame()
{
    SMCSystem::beforeEndingFrame();

	Util::scalarsVectorToOctaveFileStream(candidateChannelOrders,"candidateOrders",_f);
	Util::scalarsVectorToOctaveFileStream(iAlgorithmsPerformingChannelOrderAPPestimation,"iAlgorithmsPerformingChannelOrderAPPestimation",_f);

	channelOrderAPPsAlongTime.push_back(presentFrameChannelOrderAPPsAlongTime);
	Util::matricesVectorsVectorsVectorToOctaveFileStream(channelOrderAPPsAlongTime,"channelOrderAPPsAlongTime",_f);
}

void ChannelOrderEstimationSystem::onlyOnce()
{
	SMCSystem::onlyOnce();

	// we find out which algorithms perform channel order APP estimation
	for(uint iAlgorithm=0;iAlgorithm<_algorithms.size();iAlgorithm++)
	{
		if(_algorithms[iAlgorithm]->performsChannelOrderAPPestimation())
			// +1 is because in Octave/Matlab (where this information is supposed to be useful) there is no 0 index
			iAlgorithmsPerformingChannelOrderAPPestimation.push_back(iAlgorithm+1);
	}

	// we set the size of the results matrix for channel order APPs evolution according to the number of algorithms counted above
	presentFrameChannelOrderAPPsAlongTime = vector<vector<MatrixXd> >(iAlgorithmsPerformingChannelOrderAPPestimation.size(),vector<MatrixXd>(_SNRs.size(),MatrixXd::Zero(candidateChannelOrders.size(),_frameLength)));
}

void ChannelOrderEstimationSystem::beforeEndingAlgorithm()
{
	SMCSystem::beforeEndingAlgorithm();

	if(_algorithms[_iAlgorithm]->performsChannelOrderAPPestimation())
	{
		//...the probability of the different channel orders at each time instant is retrieved
		presentFrameChannelOrderAPPsAlongTime[iAlgorithmPerformingChannelOrderAPPestimation][_iSNR] = (dynamic_cast <UnknownChannelOrderAlgorithm *>(_algorithms[_iAlgorithm]))->getChannelOrderAPPsAlongTime_eigen();

		iAlgorithmPerformingChannelOrderAPPestimation++;
	}
}

void ChannelOrderEstimationSystem::addAlgorithms()
{
	iAlgorithmPerformingChannelOrderAPPestimation = 0;
}