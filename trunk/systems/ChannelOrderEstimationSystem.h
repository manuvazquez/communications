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
#ifndef CHANNELORDERESTIMATIONSYSTEM_H
#define CHANNELORDERESTIMATIONSYSTEM_H

#include <SMCSystem.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <ChannelOrderEstimatorSMCAlgorithm.h>

class ChannelOrderEstimationSystem : public SMCSystem
{
private:
	int iAlgorithmPerformingChannelOrderAPPestimation;
protected:
	vector<int> candidateChannelOrders;
	int iTrueChannelOrder;

	vector<tMatrix> channelOrderCoefficientsMeans;
	vector<tMatrix> channelOrderCoefficientsVariances;

	// channel order APP evolution
    vector<vector<vector<tMatrix> > > channelOrderAPPsAlongTime;
    vector<vector<tMatrix> > presentFrameChannelOrderAPPsAlongTime;
    vector<int> iAlgorithmsPerformingChannelOrderAPPestimation;

	virtual void OnlyOnce();
    virtual void BeforeEndingAlgorithm(int iAlgorithm);
	virtual void AddAlgorithms();
    virtual void BeforeEndingFrame(int iFrame);
public:
    ChannelOrderEstimationSystem();
};

#endif
