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

// #define PRINT_INFO

ChannelDependentNoise::ChannelDependentNoise(MIMOChannel *channel)
 : Noise(channel->nOutputs(),channel->length()),_matrix(StatUtil::randnMatrix(_nOutputs,_length,0.0,1.0)),_channel(channel)
{
    _stdDevs = new double[_length];
    for(int i=0;i<_length;i++)
        _stdDevs[i] = 1.0;
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

void ChannelDependentNoise::setSNR(int SNR,double alphabetVariance)
{
    int i,j,k;
    double varianceConstant = pow(10.0,((double)-SNR)/10.0)*alphabetVariance/_nOutputs;
    double stdDev,variance;

    for(j=_channel->effectiveMemory()-1;j<_length;j++)
    {
        tMatrix channelMatrix;
        channelMatrix = _channel->getTransmissionMatrix(j);
        
        variance = 0.0;
        for(i=0;i<channelMatrix.rows();i++)
            for(k=0;k<channelMatrix.cols();k++)
                variance += channelMatrix(i,k)*channelMatrix(i,k);

        variance *= varianceConstant;
        stdDev = sqrt(variance);

#ifdef PRINT_INFO
        cout << "_nOutputs = " << _nOutputs << endl;
        cout << "stdDev en el instante " << j << " = " << stdDev << endl;
#endif

        for(i=0;i<_nOutputs;i++)
        {
            // normalize by dividing by the old standard deviation and multiplying by the new one
            _matrix(i,j) = _matrix(i,j)/_stdDevs[j]*stdDev;
        }
        _stdDevs[j] = stdDev;
    }
}

double ChannelDependentNoise::stdDevAt(int n) const
{
    return _stdDevs[n];
}

tVector ChannelDependentNoise::operator[](int n) const
{
    tVector res(_nOutputs);
    for(int i=0;i<_nOutputs;i++)
        res(i) = _matrix(i,n);
    return res;
}

tMatrix ChannelDependentNoise::range(int start,int end) const
{
    tMatrix res(_nOutputs,end-start+1);
    for(int i=start;i<=end;i++)
        for(int j=0;j<_nOutputs;j++)
            res(j,i-start) = _matrix(j,i);
    return res;
}
