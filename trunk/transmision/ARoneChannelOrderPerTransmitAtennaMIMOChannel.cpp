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

ARoneChannelOrderPerTransmitAtennaMIMOChannel::ARoneChannelOrderPerTransmitAtennaMIMOChannel(int nTx, int nRx, int length, const std::vector< int >& candidateOrders, const tMatrix& channelOrderMatrixProbabilities,double mean,double variance,vector<double> ARcoefficients,double ARvariance): OneChannelOrderPerTransmitAtennaMIMOChannel(nTx, nRx, length, candidateOrders, channelOrderMatrixProbabilities)
//ARprocess constructor call
,_ARproc(StatUtil::RandnMatrix(nRx,_nChannelMatrixNotNullColumns,mean,variance),
ARcoefficients,ARvariance)
,_channelMatrices(new tMatrix[_length])
{
	tMatrix arProcessMatrix;
	tMatrix channelMatrix = LaGenMatDouble::zeros(_nRx,_nTx*_maxChannelOrder);
	int nextNotUsedColumn;
	tRange rAllRows(0,_nRx-1);

	for(int i=_maxChannelOrder-1;i<_length;i++)
	{
		arProcessMatrix = _ARproc.NextMatrix();
		nextNotUsedColumn = 0;
		for(int iAntenna=0;iAntenna<_nTx;iAntenna++)
		{
			tRange rChannelMatrixColumns((_maxChannelOrder-_antennaChannelOrder[iAntenna])*_nTx+iAntenna,_nTx*_maxChannelOrder-1,_nTx);

			#ifdef DEBUG
				cout << "El rango " << endl << rChannelMatrixColumns << endl;
			#endif

			channelMatrix(rAllRows,rChannelMatrixColumns).inject(arProcessMatrix(rAllRows,tRange(nextNotUsedColumn,nextNotUsedColumn+_antennaChannelOrder[iAntenna]-1)));

			#ifdef DEBUG
				cout << "Matriz del proceso AR" << endl << arProcessMatrix << endl;
				cout << "Matrix de canal" << endl << channelMatrix << endl;
			#endif

			nextNotUsedColumn += _antennaChannelOrder[iAntenna];

			#ifdef DEBUG
				cout << "next....: " << nextNotUsedColumn << endl;
			#endif
		}
	}
}


ARoneChannelOrderPerTransmitAtennaMIMOChannel::~ARoneChannelOrderPerTransmitAtennaMIMOChannel()
{
	delete[] _channelMatrices;
}


