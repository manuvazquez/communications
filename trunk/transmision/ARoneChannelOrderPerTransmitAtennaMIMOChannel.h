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
#ifndef ARONECHANNELORDERPERTRANSMITATENNAMIMOCHANNEL_H
#define ARONECHANNELORDERPERTRANSMITATENNAMIMOCHANNEL_H

#include <OneChannelOrderPerTransmitAtennaMIMOChannel.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <ARprocess.h>

class ARoneChannelOrderPerTransmitAtennaMIMOChannel : public OneChannelOrderPerTransmitAtennaMIMOChannel
{
private:
	inline void ConstructorProcessing();
protected:
	tMatrix* _channelMatrices;
	ARprocess _ARproc;
public:
    ARoneChannelOrderPerTransmitAtennaMIMOChannel(int nTx, int nRx, int length, const std::vector< int >& candidateOrders, const tMatrix& channelOrderMatrixProbabilities,double mean,double variance,vector<double> ARcoefficients,double ARvariance);

    ARoneChannelOrderPerTransmitAtennaMIMOChannel(int nTx, int nRx, int length,const std::vector<int> &antennasChannelOrders,double mean,double variance,vector<double> ARcoefficients,double ARvariance);

	tMatrix& operator[](int n) const { return _channelMatrices[n];};

    ~ARoneChannelOrderPerTransmitAtennaMIMOChannel();

};

#endif