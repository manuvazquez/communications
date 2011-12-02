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
#include "StillMemoryMIMOChannel.h"

#include <assert.h>

StillMemoryMIMOChannel::StillMemoryMIMOChannel(uint nInputs, uint nOutputs, uint memory, uint length): MIMOChannel(nInputs, nOutputs, length),_memory(memory),_nInputsnOutputsMemory(_nInputs*_nOutputs*_memory),_nInputsMemory(_nInputs*_memory)
{
}

void StillMemoryMIMOChannel::setSubchannelOrders(std::vector<uint> subchannelOrders)
{
  if(subchannelOrders.size()!=static_cast<uint>(_nOutputs))
	throw RuntimeException("StillMemoryMIMOChannel::setSubchannelOrders: number of channel orders is not equal to the number of outputs.");
  
  MatrixXd matrix;
  uint maxChannelOrderFound=0;

  // we need to check that...
  for(uint i=0;i<subchannelOrders.size();i++)
  {
	// ...the overall channel order is present...
	maxChannelOrderFound += (subchannelOrders[i]==static_cast<uint>(_memory));
	
	//...and that all the channel orders passed are smaller than the overall channel order
	if(subchannelOrders[i]>static_cast<uint>(_memory))
	  throw RuntimeException("StillMemoryMIMOChannel::setSubchannelOrders: one subchannel order or more is bigger than the overall channel order.");
  }
  
//   if(!maxChannelOrderFound)
// 	throw RuntimeException("StillMemoryMIMOChannel::setSubchannelOrders: none of the subchannel orders is equal to the maximum.");
  
  // FIXME: this implementation is inefficient when the channel is "TimeInvariantChannel"
  
  assert(_memory>0);
  for(uint iMat=_memory-1;iMat<_length;iMat++)
  {
	matrix = at(iMat);
	for(uint iRow=0;iRow<static_cast<uint>(_nOutputs);iRow++)
	  for(uint iCol=0;iCol<static_cast<uint>(matrix.cols())-_nInputs*subchannelOrders[iRow];iCol++)
		matrix(iRow,iCol) = 0.0;
	set(iMat,matrix);
  }
}
