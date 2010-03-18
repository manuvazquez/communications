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
#include "ARchannel.h"

using namespace std;

ARchannel::ARchannel(int nInputs, int nOutputs, int memory, int length, ARprocess ARproc): StillMemoryMIMOChannel(nInputs, nOutputs, memory, length),_channelMatrices(length),_ARproc(ARproc)
{
	if(ARproc.rows()!=nOutputs || ARproc.cols()!=(nInputs*memory))
		throw RuntimeException("ARchannel::ARchannel: the passed AR process is not compatible with the dimensions of the channel.");

	//initialization
	for(int i=_memory-1;i<_length;i++)
			_channelMatrices[i] = _ARproc.nextMatrix();
}
