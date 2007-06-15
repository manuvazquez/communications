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
#include "ARoneChannelOrderPerTransmitAtennaMIMOChannel.h"

// #define DEBUG

ARoneChannelOrderPerTransmitAtennaMIMOChannel::ARoneChannelOrderPerTransmitAtennaMIMOChannel(int nTx, int nRx, int length,const std::vector<int> &antennasChannelOrders,double mean,double variance,vector<double> ARcoefficients,double ARvariance): OneChannelOrderPerTransmitAtennaMIMOChannel(nTx, nRx, length,antennasChannelOrders)
//ARprocess constructor call
,_channelMatrices(new tMatrix[_length]),_ARproc(StatUtil::RandnMatrix(nRx,_nChannelMatrixNotNullColumns,mean,variance),
ARcoefficients,ARvariance)
{
	#ifdef DEBUG
		cout << "Antes de empezar el bucle para crear las matrice para los instantes de tiempo." << endl;
	#endif

	for(int i=_memory-1;i<_length;i++)
	{
		tMatrix channelMatrix(_nRx,_nTx*_memory);
		tMatrix arProcessMatrix = _ARproc.NextMatrix();

		OneChannelOrderPerTransmitAtennaMIMOChannel::WithoutZerosMatrixToWithZerosMatrix(arProcessMatrix,_nTx,_antennasChannelOrders,channelMatrix);

		_channelMatrices[i] = channelMatrix;
	}

	#ifdef DEBUG
		cout << "Despues del bucle para crear las matrice para los instantes de tiempo." << endl;
	#endif
}

ARoneChannelOrderPerTransmitAtennaMIMOChannel::~ARoneChannelOrderPerTransmitAtennaMIMOChannel()
{
	delete[] _channelMatrices;
}
