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
#ifndef STILLMEMORYMIMOCHANNEL_H
#define STILLMEMORYMIMOCHANNEL_H

#include <MIMOChannel.h>

/**
	@author Manu <manu@rustneversleeps>
*/
class StillMemoryMIMOChannel : public MIMOChannel
{
protected:
	uint _memory,_nInputsnOutputsMemory,_nInputsMemory;
public:
    StillMemoryMIMOChannel(uint nInputs, uint nOutputs, uint memory,uint length);

	uint memory() const {return _memory;};
	uint memory(uint n) const {return _memory;}
	uint effectiveMemory() const {return _memory;}
	uint nInputsnOutputsmemory() const {return _nInputsnOutputsMemory;};
	uint nInputsmemory() const {return _nInputsMemory;};
	
	//! it shapes the channel matrices of the channel so that the outputs have the specified channel orders (introducing zeros whenever necessary)
	/*!
		\param channelOrders a vector of uint with the channel orders of the different outputs
	*/
	void setSubchannelOrders(std::vector<uint> subchannelOrders);
};

#endif
