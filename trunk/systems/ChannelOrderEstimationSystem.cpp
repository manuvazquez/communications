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

ChannelOrderEstimationSystem::ChannelOrderEstimationSystem()
 : SMCSystem()
{
	candidateChannelOrders.push_back(2);candidateChannelOrders.push_back(3);candidateChannelOrders.push_back(4);
	candidateChannelOrders.push_back(5);candidateChannelOrders.push_back(6);candidateChannelOrders.push_back(7);

	channelOrderCoefficientsMeans.resize(candidateChannelOrders.size());
	channelOrderCoefficientsVariances.resize(candidateChannelOrders.size());

	iTrueChannelOrder = -1;
	for(uint iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		channelOrderCoefficientsMeans[iChannelOrder] = LaGenMatDouble::zeros(L,N*candidateChannelOrders[iChannelOrder]);
		channelOrderCoefficientsVariances[iChannelOrder] = LaGenMatDouble::ones(L,N*candidateChannelOrders[iChannelOrder]);
		if(candidateChannelOrders[iChannelOrder] == m)
			iTrueChannelOrder = iChannelOrder;
	}

	if(iTrueChannelOrder==-1)
		throw RuntimeException("ChannelOrderEstimationSystem::ChannelOrderEstimationSystem: the memory of the channel is not one of the possible candidates.");
}

void ChannelOrderEstimationSystem::BeforeEndingFrame(int iFrame)
{
	Util::ScalarsVectorToOctaveFileStream(candidateChannelOrders,"candidateOrders",f);
    SMCSystem::BeforeEndingFrame(iFrame);
}

