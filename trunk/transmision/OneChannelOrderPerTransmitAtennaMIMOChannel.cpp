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
#include "OneChannelOrderPerTransmitAtennaMIMOChannel.h"

#define DEBUG

OneChannelOrderPerTransmitAtennaMIMOChannel::OneChannelOrderPerTransmitAtennaMIMOChannel(int nTx, int nRx, int length,const std::vector<int>& candidateOrders,const tMatrix& channelOrderMatrixProbabilities): MIMOChannel(nTx, nRx, length)
,_channelOrderMatrixProbabilities(channelOrderMatrixProbabilities),_candidateOrders(candidateOrders),_antennaChannelOrder(new int[_nTx])
{
	_maxChannelOrder = -1;
	_nChannelMatrixNotNullColumns = 0;
	for(int i=0;i<_nTx;i++)
	{
		_antennaChannelOrder[i] = _candidateOrders[StatUtil::Discrete_rnd(channelOrderMatrixProbabilities.row(i))];

		if(_antennaChannelOrder[i] > _maxChannelOrder)
			_maxChannelOrder = _antennaChannelOrder[i];

		_nChannelMatrixNotNullColumns += _antennaChannelOrder[i];
	}

	#ifdef DEBUG
		for(int i=0;i<_nTx;i++)
			cout << _antennaChannelOrder[i] << endl;

		cout << "El maximo: " << _maxChannelOrder << " y el nº total de columnas: " << _nChannelMatrixNotNullColumns << endl;
	#endif
}


OneChannelOrderPerTransmitAtennaMIMOChannel::~OneChannelOrderPerTransmitAtennaMIMOChannel()
{
	delete[] _antennaChannelOrder;
}
