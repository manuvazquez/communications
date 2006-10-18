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
#ifndef SPRAWLINGMEMORYARMIMOCHANNEL_H
#define SPRAWLINGMEMORYARMIMOCHANNEL_H

#include <SprawlingMemoryMIMOChannel.h>
#include <ARprocess.h>

/**
	@author Manu <manu@rustneversleeps>
*/
class SprawlingMemoryARMIMOChannel : public SprawlingMemoryMIMOChannel
{
protected:
	tMatrix* _channelMatrices;
	int *_channelOrders;
	ARprocess _ARprocess;
public:
    SprawlingMemoryARMIMOChannel(int nTx, int nRx, int length, std::vector< int > candidateOrders, tMatrix transitionProbabilitiesMatrix, int initialChannelOrderIndex,double mean,double variance,std::vector<double> ARcoefficients,double ARvariance,Random randomGenerator =  Random(0));

    ~SprawlingMemoryARMIMOChannel();

	int Memory(int n) const {return _channelOrders[n];}
	tMatrix& operator[](int n) const { return _channelMatrices[n];};
};

#endif
