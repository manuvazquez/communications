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

ChannelDependentNoise::ChannelDependentNoise(MIMOChannel &channel)
 : Noise(channel.Nr(),channel.Length()),channel(channel)
{
	stdDevs = new double[length];
	for(int i=0;i<length;i++)
		stdDevs[i] = 1;
}


ChannelDependentNoise::~ChannelDependentNoise()
{
}

// 		public void EstablecerSNR(int SNR,double varianzaAlfabeto)
// 		{
// 			int i,j,memoriaCanal;
// 			double cteVarianza = Math.Pow(10d,((double)-SNR)/10d);
// 			double desviacionTipica,varianza;
//
// 			memoriaCanal = canal.Memoria;
//
// 			for(j=memoriaCanal-1;j<nInstantesTiempo;j++)
// 			{
// 				varianza = cteVarianza*((canal[j].Transpose()*canal[j]).Diagonal()).Sum()*//
// 						   varianzaAlfabeto/nAntenasReceptoras;
// 				desviacionTipica = Math.Sqrt(varianza);
// 				for(i=0;i<nAntenasReceptoras;i++)
// 				{
// 					// normaliza dividiendo por la desviacion tipica anterior, y multiplica
// 					// por la nueva
// 					ruido[i,j] = ruido[i,j]/desviacionesTipicas[j]*desviacionTipica;
// 				}
// 				desviacionesTipicas[j] = desviacionTipica;
// 			}
// 		}

void ChannelDependentNoise::SetSNR(int SNR,double alphabetVariance)
{
	int i,j,channelMemory;
	double varianceConstant = pow(10.0,((double)-SNR)/10.0);
	double stdDev,variance;

	channelMemory = channel.Memory();

	tMatrix channelTranspChannel(channel.NtMemory(),channel.NtMemory());
	for(j=channelMemory-1;j<length;j++)
	{
		//channelTranspChannel = varianceConstant*alphabetVariance/Nr*channel[j]'*channel[j];
		Blas_Mat_Trans_Mat_Mult(channel[j],channel[j],channelTranspChannel,varianceConstant*alphabetVariance/nRx);

		variance = channelTranspChannel.trace();
		stdDev = sqrt(variance);

		for(i=0;i<nRx;i++)
			// normalize by dividing by the old standard deviation and multiplying by the new one
			matrix(i,j) = matrix(i,j)/stdDevs[j]*stdDev;
		stdDevs[j] = stdDev;
	}
}

double ChannelDependentNoise::StdDevAt(int n)
{
	return stdDevs[n];
}

tVector ChannelDependentNoise::operator[](int n)
{
	tVector res(nRx);
	for(int i=0;i<nRx;i++)
		res(i) = matrix(i,n);
	return res;
}
