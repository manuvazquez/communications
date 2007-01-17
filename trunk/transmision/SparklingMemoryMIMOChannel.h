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
#ifndef SPARKLINGMEMORYMIMOCHANNEL_H
#define SPARKLINGMEMORYMIMOCHANNEL_H

#include <MIMOChannel.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <Util.h>

class SparklingMemoryMIMOChannel : public MIMOChannel
{
protected:
	std::vector<int> _candidateOrders;
	int _maxOrder,_initialChannelOrderIndex;
public:
    SparklingMemoryMIMOChannel(int nTx,int nRx,int length,std::vector<int> candidateOrders,int initialChannelOrderIndex);

    int EffectiveMemory() const {return _maxOrder;}
};

#endif
