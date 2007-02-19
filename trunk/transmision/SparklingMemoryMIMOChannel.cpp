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
#include "SparklingMemoryMIMOChannel.h"

using namespace std;

SparklingMemoryMIMOChannel::SparklingMemoryMIMOChannel(int nTx, int nRx, int length, vector<int> candidateOrders,int initialChannelOrderIndex): MIMOChannel(nTx, nRx, length),_candidateOrders(candidateOrders),_initialChannelOrderIndex(initialChannelOrderIndex)
{
	// it obtains the maximum allowed order for the channel
	int iMaxOrder = 0;
	_maxOrder = _candidateOrders[0];
	for(uint i=1;i<_candidateOrders.size();i++)
		if(_candidateOrders[i]>_maxOrder)
		{
			iMaxOrder = i;
			_maxOrder = _candidateOrders[i];
		}
}
