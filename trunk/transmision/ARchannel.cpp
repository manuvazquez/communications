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

ARchannel::ARchannel(int nTx, int nRx, int memory, int length,double mean,double variance,vector<double> ARcoefficients,double ARvariance,Random randomGenerator): MIMOChannel(nTx, nRx, memory, length),
//ARprocess constructor call
_ARproc(StatUtil::RandnMatrix(nRx,nTx*memory,mean,variance,randomGenerator),
ARcoefficients,ARvariance,randomGenerator)
{
	_channelMatrices = new tMatrix[length];

	//initialization
	for(int i=_memory-1;i<_length;i++)
			_channelMatrices[i] = _ARproc.NextMatrix();
}

ARchannel::ARchannel(const ARchannel &archannel):MIMOChannel(archannel),_ARproc(archannel._ARproc),_channelMatrices(new tMatrix[_length])
{
	for(int i=_memory-1;i<_length;i++)
			_channelMatrices[i] = archannel._channelMatrices[i];
}

ARchannel::~ ARchannel()
{
	delete[] _channelMatrices;
}

vector<tMatrix> ARchannel::Range(int a,int b)
{
	int nMatrices = b - a + 1;

	if(nMatrices<1)
		throw RuntimeException("Selected range of time is invalid.");
	vector<tMatrix> res(nMatrices);

	for(int i=0;i<nMatrices;i++)
	{
		res[i] = _channelMatrices[a+i];
// 		cout << "en ARChannel" << endl << _channelMatrices[i] << endl;
	}

	return res;
}
