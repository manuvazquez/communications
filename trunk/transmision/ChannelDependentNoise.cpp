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
#include "ChannelDependentNoise.h"

// #define DEBUG

ChannelDependentNoise::ChannelDependentNoise(MIMOChannel *channel)
 : Noise(channel->Nr(),channel->Length()),_channel(channel)
{
	_stdDevs = new double[_length];
	for(int i=0;i<_length;i++)
		_stdDevs[i] = 1;
}

ChannelDependentNoise::ChannelDependentNoise(const ChannelDependentNoise &channelDependentNoise):Noise(channelDependentNoise),_channel(channelDependentNoise._channel),_stdDevs(new double[_length])
{
	for(int i=0;i<_length;i++)
		_stdDevs[i] = channelDependentNoise._stdDevs[i];
}

ChannelDependentNoise::~ChannelDependentNoise()
{
	delete[] _stdDevs;
}

void ChannelDependentNoise::SetSNR(int SNR,double alphabetVariance)
{
	int i,j,memory;
	double varianceConstant = pow(10.0,((double)-SNR)/10.0)*alphabetVariance/_nRx;
	double stdDev,variance;

	for(j=_channel->EffectiveMemory()-1;j<_length;j++)
	{
		memory = _channel->Memory(j);

		tMatrix channelTranspChannel(_channel->NtMemory(j),_channel->NtMemory(j));

		//channelTranspChannel = varianceConstant*alphabetVariance/Nr*_channel[j]'*_channel[j];
		Blas_Mat_Trans_Mat_Mult((*_channel)[j],(*_channel)[j],channelTranspChannel,varianceConstant);

		variance = channelTranspChannel.trace();
		stdDev = sqrt(variance);

		for(i=0;i<_nRx;i++)
		{
			// normalize by dividing by the old standard deviation and multiplying by the new one
			_matrix(i,j) = _matrix(i,j)/_stdDevs[j]*stdDev;
		}
		_stdDevs[j] = stdDev;
	}
}

double ChannelDependentNoise::StdDevAt(int n)
{
	return _stdDevs[n];
}

tVector ChannelDependentNoise::operator[](int n)
{
	tVector res(_nRx);
	for(int i=0;i<_nRx;i++)
		res(i) = _matrix(i,n);
	return res;
}

tMatrix ChannelDependentNoise::Range(int start,int end) const
{
	tMatrix res(_nRx,end-start+1);
	for(int i=start;i<=end;i++)
		for(int j=0;j<_nRx;j++)
			res(j,i-start) = _matrix(j,i);
	return res;
}
