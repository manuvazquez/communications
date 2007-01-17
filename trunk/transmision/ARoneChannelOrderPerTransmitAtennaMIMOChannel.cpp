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

void ARoneChannelOrderPerTransmitAtennaMIMOChannel::ConstructorProcessing()
{
	for(int i=_maxChannelOrder-1;i<_length;i++)
	{
		tMatrix channelMatrix(_nRx,_nTx*_maxChannelOrder);
		tMatrix arProcessMatrix = _ARproc.NextMatrix();

		OneChannelOrderPerTransmitAtennaMIMOChannel::WithoutZerosMatrixToWithZerosMatrix(arProcessMatrix,_nTx,_antennasChannelOrders,channelMatrix);

		_channelMatrices[i] = channelMatrix;
	}
}

ARoneChannelOrderPerTransmitAtennaMIMOChannel::ARoneChannelOrderPerTransmitAtennaMIMOChannel(int nTx, int nRx, int length, const std::vector< int >& candidateOrders, const tMatrix& channelOrderMatrixProbabilities,double mean,double variance,vector<double> ARcoefficients,double ARvariance): OneChannelOrderPerTransmitAtennaMIMOChannel(nTx, nRx, length, candidateOrders, channelOrderMatrixProbabilities)
//ARprocess constructor call
,_ARproc(StatUtil::RandnMatrix(nRx,_nChannelMatrixNotNullColumns,mean,variance),
ARcoefficients,ARvariance)
,_channelMatrices(new tMatrix[_length])
{
// 	for(int i=_maxChannelOrder-1;i<_length;i++)
// 	{
// 		tMatrix channelMatrix(_nRx,_nTx*_maxChannelOrder);
// 		tMatrix arProcessMatrix = _ARproc.NextMatrix();
//
// 		#ifdef DEBUG
// 			tMatrix matrizRecuperada(arProcessMatrix);
// 			matrizRecuperada = 0.0;
// 			cout << "matriz que se va a intercalar con 0s" << endl << arProcessMatrix << endl;
// 			cout << "Matriz sin ceros antes" << endl << matrizRecuperada << endl;
// 		#endif
//
// 		OneChannelOrderPerTransmitAtennaMIMOChannel::WithoutZerosMatrixToWithZerosMatrix(arProcessMatrix,_nTx,_antennasChannelOrders,channelMatrix);
//
// 		#ifdef DEBUG
// 			cout << "Matriz intercalada" << endl << channelMatrix << endl;
// 		#endif
//
// 		#ifdef DEBUG
// 			OneChannelOrderPerTransmitAtennaMIMOChannel::WithZerosMatrixToWithoutZerosMatrix(channelMatrix,_nTx,antennasChannelOrders,matrizRecuperada);
// 			cout << "Matriz sin ceros despues de recuperada" << endl << matrizRecuperada << endl;
// 			getchar();
// 		#endif
//
// 		_channelMatrices[i] = channelMatrix;
// 	}

	ConstructorProcessing();
}

ARoneChannelOrderPerTransmitAtennaMIMOChannel::ARoneChannelOrderPerTransmitAtennaMIMOChannel(int nTx, int nRx, int length,const std::vector<int> &antennasChannelOrders,double mean,double variance,vector<double> ARcoefficients,double ARvariance): OneChannelOrderPerTransmitAtennaMIMOChannel(nTx, nRx, length,antennasChannelOrders)
//ARprocess constructor call
,_ARproc(StatUtil::RandnMatrix(nRx,_nChannelMatrixNotNullColumns,mean,variance),
ARcoefficients,ARvariance)
,_channelMatrices(new tMatrix[_length])
{
	ConstructorProcessing();
}

ARoneChannelOrderPerTransmitAtennaMIMOChannel::~ARoneChannelOrderPerTransmitAtennaMIMOChannel()
{
	delete[] _channelMatrices;
}
