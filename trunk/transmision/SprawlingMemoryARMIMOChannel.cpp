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
#include "SprawlingMemoryARMIMOChannel.h"

SprawlingMemoryARMIMOChannel::SprawlingMemoryARMIMOChannel(int nTx, int nRx, int length, std::vector< int > candidateOrders, tMatrix transitionProbabilitiesMatrix, int initialChannelOrderIndex,double mean,double variance,vector<double> ARcoefficients,double ARvariance,Random randomGenerator): SprawlingMemoryMIMOChannel(nTx, nRx, length, candidateOrders, transitionProbabilitiesMatrix, initialChannelOrderIndex),_ARprocess(StatUtil::RandnMatrix(nRx,nTx*_maxOrder,mean,variance),
ARcoefficients,ARvariance,randomGenerator)
{
	int presentChannelOrder,iPresentChannelOrder;
	tRange rAllObservationsRows(0,_nRx-1);

	_channelMatrices = new tMatrix[length];
	_channelOrders = new int[length];

	// the first channel matrix has memory "initialChannelOrderIndex"
	iPresentChannelOrder = initialChannelOrderIndex;
	presentChannelOrder = _candidateOrders[_initialChannelOrderIndex];

	for(int i=_maxOrder-1;i<_length;i++)
	{
		_channelOrders[i] = iPresentChannelOrder;

		_channelMatrices[i] = _ARprocess.NextMatrix()(rAllObservationsRows,tRange(0,_nTx*_candidateOrders[iPresentChannelOrder]-1));

		iPresentChannelOrder = StatUtil::Discrete_rnd(1,transitionProbabilitiesMatrix.row(iPresentChannelOrder))[0];
	}

	// the first channel order is modified in the loop but it must be "initialChannelOrderIndex"
	_channelOrders[_maxOrder-1] = initialChannelOrderIndex;
}

SprawlingMemoryARMIMOChannel::SprawlingMemoryARMIMOChannel(const SprawlingMemoryARMIMOChannel &channel):SprawlingMemoryMIMOChannel(channel),_ARprocess(channel._ARprocess),_channelMatrices(new tMatrix[_length]),_channelOrders(new int[_length])
{
	for(int i=_maxOrder-1;i<_length;i++)
	{
		_channelOrders[i] = channel._channelOrders[i];
		_channelMatrices[i] = channel._channelMatrices[i];
	}
}

SprawlingMemoryARMIMOChannel::~SprawlingMemoryARMIMOChannel()
{
	delete[] _channelMatrices;
	delete[] _channelOrders;
}

