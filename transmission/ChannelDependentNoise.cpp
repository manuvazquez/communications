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

#include <assert.h>
// #define PRINT_INFO

ChannelDependentNoise::ChannelDependentNoise(double alphabetVariance,MIMOChannel *channel)
 : Noise(channel->nOutputs(),channel->length()),_matrix(StatUtil::randnMatrix(_nOutputs,_length,0.0,1.0)),_channel(channel),_stdDevs(_length),_alphabetVariance(alphabetVariance)
{
    for(uint i=0;i<_length;i++)
        _stdDevs[i] = 1.0;
}

void ChannelDependentNoise::setSNR(int SNR)
{
    uint i,j;
    double SNRdependentVarianceFactor = pow(10.0,((double)-SNR)/10.0)*_alphabetVariance/_nOutputs;
    double stdDev;

	assert(_channel->effectiveMemory()>0);
    for(j=_channel->effectiveMemory()-1;j<_length;j++)
    {
        MatrixXd channelMatrix = _channel->getTransmissionMatrix(j);
		
		stdDev = computeStd(SNRdependentVarianceFactor,channelMatrix);

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

VectorXd ChannelDependentNoise::at(uint n) const
{
    return _matrix.col(n);
}

MatrixXd ChannelDependentNoise::range(int start,int end) const
{
    return _matrix.block(0,start,_nOutputs,end-start+1);
}
