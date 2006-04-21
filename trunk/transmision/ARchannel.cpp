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

ARchannel::ARchannel(int nTx, int nRx, int memory, int length,double mean,double variance,vector<double> ARcoefficients,double ARvariance,Random &randomGenerator): MIMOChannel(nTx, nRx, memory, length),
//ARprocess constructor call
// ARproc(*(new tMatrix(randomGenerator.randnArray(nTx*nRx*memory,mean,variance),nRx,nTx*memory)),
// ARcoefficients,ARvariance)
ARproc(StatUtil::RandnMatrix(nRx,nTx*memory,mean,variance,randomGenerator),
ARcoefficients,ARvariance)

{
	channelMatrices = new tMatrix[length];

	//initialization
	for(int i=memory-1;i<length;i++)
		channelMatrices[i] = ARproc.NextMatrix();
}

ARchannel::~ ARchannel()
{
	delete[] channelMatrices;
}
