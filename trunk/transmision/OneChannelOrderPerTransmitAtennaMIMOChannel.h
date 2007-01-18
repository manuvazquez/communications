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
#ifndef ONECHANNELORDERPERTRANSMITATENNAMIMOCHANNEL_H
#define ONECHANNELORDERPERTRANSMITATENNAMIMOCHANNEL_H

#include <StillMemoryMIMOChannel.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <StatUtil.h>

class OneChannelOrderPerTransmitAtennaMIMOChannel : public StillMemoryMIMOChannel
{
protected:
	std::vector<int> _antennasChannelOrders;
	int _nChannelMatrixNotNullColumns;
public:
    OneChannelOrderPerTransmitAtennaMIMOChannel(int nTx, int nRx, int length,const std::vector<int> &antennasChannelOrders);

    vector<int> GetAntennasChannelOrders() {return _antennasChannelOrders;}

	static void WithoutZerosMatrixToWithZerosMatrix(const tMatrix &withoutZerosMatrix,int N,const vector<int> &antennasChannelOrders,tMatrix &withZerosMatrix);

	static tMatrix WithZerosMatrixToWithoutZerosMatrix(const tMatrix &withZerosMatrix,int N,const vector<int> &antennasChannelOrders);

	static void CompleteSymbolsVectorToOnlyInvolvedSymbolsVector(const tVector &withZerosSymbolsVector,int N,const vector<int> &antennasChannelOrders,tVector &withoutZerosSymbolsVector);
};

#endif
