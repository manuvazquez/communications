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

#include<algorithm>

ChannelOrderEstimationSystem::ChannelOrderEstimationSystem()
 : SMCSystem()
{
	_candidateChannelOrders.push_back(1);
	_candidateChannelOrders.push_back(2);
	_candidateChannelOrders.push_back(3);
	_candidateChannelOrders.push_back(4);
	
// 	_candidateChannelOrders.push_back(5);
// 	_candidateChannelOrders.push_back(6);
// 	_candidateChannelOrders.push_back(7);

	_channelOrderCoefficientsMeans.resize(_candidateChannelOrders.size());
	_channelOrderCoefficientsVariances.resize(_candidateChannelOrders.size());

	_iMaxChannelOrder = std::max_element(_candidateChannelOrders.begin(),_candidateChannelOrders.end())-_candidateChannelOrders.begin();
	
	std::vector<uint>::iterator trueChanelOrderIterator = std::find(_candidateChannelOrders.begin(),_candidateChannelOrders.end(),_m);
	
	if(trueChanelOrderIterator== _candidateChannelOrders.end())
		throw RuntimeException("ChannelOrderEstimationSystem::ChannelOrderEstimationSystem: the memory of the channel is not one of the possible candidates.");
	
	_iTrueChannelOrder = trueChanelOrderIterator - _candidateChannelOrders.begin();
	
	for(uint iChannelOrder=0;iChannelOrder<_candidateChannelOrders.size();iChannelOrder++)
	{
		_channelOrderCoefficientsMeans[iChannelOrder] = MatrixXd::Zero(_L,_N*_candidateChannelOrders[iChannelOrder]);
		_channelOrderCoefficientsVariances[iChannelOrder] = MatrixXd::Ones(_L,_N*_candidateChannelOrders[iChannelOrder]);
	}

	// channel order APP evolution
    _channelOrderAPPsAlongTime.reserve(_nFrames);
}

void ChannelOrderEstimationSystem::onlyOnce()
{
	SMCSystem::onlyOnce();

	// we find out which algorithms perform channel order APP estimation
	for(uint iAlgorithm=0;iAlgorithm<_algorithms.size();iAlgorithm++)
	{
		if(_algorithms[iAlgorithm]->estimatesOneSingleChannelOrder())
			// +1 is because in Octave/Matlab (where this information is supposed to be useful) there is no 0 index
			_iAlgorithmsPerformingChannelOrderAPPestimation.push_back(iAlgorithm+1);
		
		if(_algorithms[iAlgorithm]->estimatesOneChannelOrderPerOutput())
			// +1 is because in Octave/Matlab (where this information is supposed to be useful) there is no 0 index
			_iAlgorithmsPerformingOneChannelOrderPerOutputAPPestimation.push_back(iAlgorithm+1);
	}

	// we set the size of the results matrix for channel order APPs evolution according to the number of algorithms counted above
	_presentFrameChannelOrderAPPsAlongTime = vector<vector<MatrixXd> >(_iAlgorithmsPerformingChannelOrderAPPestimation.size(),vector<MatrixXd>(_SNRs.size(),MatrixXd::Zero(_candidateChannelOrders.size(),_frameLength)));
	_presentFrameOneChannelOrderPerOutputAPPsAlongTime = vector<vector<vector<MatrixXd> > >(_L,
														  vector<vector<MatrixXd> >(_iAlgorithmsPerformingOneChannelOrderPerOutputAPPestimation.size(),
																vector<MatrixXd>(_SNRs.size(),MatrixXd::Zero(_candidateChannelOrders.size(),_frameLength))));
}

void ChannelOrderEstimationSystem::beforeEndingAlgorithm()
{
	SMCSystem::beforeEndingAlgorithm();

	if(_algorithms[_iAlgorithm]->estimatesOneSingleChannelOrder())
	{
		//...the probability of the different channel orders at each time instant is retrieved
		_presentFrameChannelOrderAPPsAlongTime[_iAlgorithmPerformingChannelOrderAPPestimation][_iSNR] = (dynamic_cast <UnknownChannelOrderAlgorithm *>(_algorithms[_iAlgorithm]))->getComputedChannelOrderAPPs();

		_iAlgorithmPerformingChannelOrderAPPestimation++;
	}
	else
	  if(_algorithms[_iAlgorithm]->estimatesOneChannelOrderPerOutput())
	  {
		for(uint iOutput=0;iOutput<_L;iOutput++)
			_presentFrameOneChannelOrderPerOutputAPPsAlongTime[iOutput][_iAlgorithmPerformingOneChannelOrderPerOutputAPPestimation][_iSNR] = (dynamic_cast <UnknownChannelOrderAlgorithm *>(_algorithms[_iAlgorithm]))->getComputedChannelOrderAPPs(iOutput);
		
		_iAlgorithmPerformingOneChannelOrderPerOutputAPPestimation++;
	  }
}

void ChannelOrderEstimationSystem::addAlgorithms()
{
	_iAlgorithmPerformingChannelOrderAPPestimation = 0;
	_iAlgorithmPerformingOneChannelOrderPerOutputAPPestimation = 0;
}

void ChannelOrderEstimationSystem::storeFrameResults()
{
    SMCSystem::storeFrameResults();

	_channelOrderAPPsAlongTime.push_back(_presentFrameChannelOrderAPPsAlongTime);
	_oneChannelOrderPerOutputAPPsAlongTime.push_back(_presentFrameOneChannelOrderPerOutputAPPsAlongTime);
}

void ChannelOrderEstimationSystem::saveFrameResults()
{
    SMCSystem::saveFrameResults();

	Octave::toOctaveFileStream(_candidateChannelOrders,"candidateOrders",_f);
	Octave::toOctaveFileStream(_iAlgorithmsPerformingChannelOrderAPPestimation,"iAlgorithmsPerformingChannelOrderAPPestimation",_f);
	Octave::toOctaveFileStream(_iAlgorithmsPerformingOneChannelOrderPerOutputAPPestimation,"iAlgorithmsPerformingOneChannelOrderPerOutputAPPestimation",_f);
	Octave::eigenToOctaveFileStream(_channelOrderAPPsAlongTime,"channelOrderAPPsAlongTime",_f);	
	Octave::eigenToOctaveFileStream(_oneChannelOrderPerOutputAPPsAlongTime,"oneChannelOrderPerOutputAPPsAlongTime",_f);
}