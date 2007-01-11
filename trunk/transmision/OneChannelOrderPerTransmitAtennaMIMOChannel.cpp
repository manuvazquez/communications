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

// #define DEBUG

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

		cout << "El maximo: " << _maxChannelOrder << " y el n� total de columnas: " << _nChannelMatrixNotNullColumns << endl;
	#endif
}


OneChannelOrderPerTransmitAtennaMIMOChannel::~OneChannelOrderPerTransmitAtennaMIMOChannel()
{
	delete[] _antennaChannelOrder;
}

void OneChannelOrderPerTransmitAtennaMIMOChannel::WithoutZerosMatrixToWithZerosMatrix(const tMatrix &withoutZerosMatrix,int N,const vector<int> &antennasChannelOrders,tMatrix &withZerosMatrix)
{
	if(withoutZerosMatrix.rows()!=withZerosMatrix.rows())
		throw RuntimeException("ARoneChannelOrderPerTransmitAtennaMIMOChannel::WithoutZerosMatrixToWithZerosMatrix: matrix dimensions are not coherent");

	int maxChannelOrder = antennasChannelOrders[Util::Max(antennasChannelOrders)];

	#ifdef DEBUG
		cout << "las columnas de la matriz destino: " << withZerosMatrix.cols() << endl;
		cout << "maxChannelOrder: " << maxChannelOrder << " N: " << N << endl;
	#endif

	if(withZerosMatrix.cols()!=(N*maxChannelOrder))
		throw RuntimeException("ARoneChannelOrderPerTransmitAtennaMIMOChannel::WithoutZerosMatrixToWithZerosMatrix: resultant matrix has not enough columns");

	tRange rAllRows(0,withoutZerosMatrix.rows()-1);

	withZerosMatrix = 0.0;

	int nextNotUsedColumn = 0;
	for(int iAntenna=0;iAntenna<N;iAntenna++)
	{
		tRange rDestColumns((maxChannelOrder-antennasChannelOrders[iAntenna])*N+iAntenna,N*maxChannelOrder-1,N);

		#ifdef DEBUG
			cout << "El rango " << endl << rDestColumns << endl;
		#endif

		withZerosMatrix(rAllRows,rDestColumns).inject(withoutZerosMatrix(rAllRows,tRange(nextNotUsedColumn,nextNotUsedColumn+antennasChannelOrders[iAntenna]-1)));

		#ifdef DEBUG
			cout << "Matriz del proceso AR" << endl << withoutZerosMatrix << endl;
			cout << "Matrix de canal" << endl << withZerosMatrix << endl;
		#endif

		nextNotUsedColumn += antennasChannelOrders[iAntenna];

		#ifdef DEBUG
			cout << "next....: " << nextNotUsedColumn << endl;
			getchar();
		#endif
	}
}

void OneChannelOrderPerTransmitAtennaMIMOChannel::WithZerosMatrixToWithoutZerosMatrix(const tMatrix &withZerosMatrix,int N,const vector<int> &antennasChannelOrders,tMatrix &withoutZerosMatrix)
{
	if(withoutZerosMatrix.rows()!=withZerosMatrix.rows())
		throw RuntimeException("ARoneChannelOrderPerTransmitAtennaMIMOChannel::WithZerosMatrixToWithoutZerosMatrix: matrix dimensions are not coherent");

	int maxChannelOrder = antennasChannelOrders[Util::Max(antennasChannelOrders)];
	int channelOrdersSum = Util::Sum(antennasChannelOrders);

	#ifdef DEBUG
		cout << "las columnas de la matriz destino: " << withZerosMatrix.cols() << endl;
		cout << "maxChannelOrder: " << maxChannelOrder << " N: " << N << endl;
	#endif

	if(withoutZerosMatrix.cols()!= channelOrdersSum)
		throw RuntimeException("ARoneChannelOrderPerTransmitAtennaMIMOChannel::WithZerosMatrixToWithoutZerosMatrix: resultant matrix has not enough columns");

	tRange rAllRows(0,withoutZerosMatrix.rows()-1);

	withoutZerosMatrix = 0.0;

	int nextNotUsedColumn = 0;
	for(int iAntenna=0;iAntenna<N;iAntenna++)
	{
		tRange rSourceColumns((maxChannelOrder-antennasChannelOrders[iAntenna])*N+iAntenna,N*maxChannelOrder-1,N);

		#ifdef DEBUG
			cout << "El rango " << endl << rSourceColumns << endl;
			cout << withoutZerosMatrix(rAllRows,tRange(nextNotUsedColumn,nextNotUsedColumn+antennasChannelOrders[iAntenna]-1)) << endl;
		#endif

		withoutZerosMatrix(rAllRows,tRange(nextNotUsedColumn,nextNotUsedColumn+antennasChannelOrders[iAntenna]-1)).inject(withZerosMatrix(rAllRows,rSourceColumns));

		#ifdef DEBUG
			cout << "Matriz del proceso AR" << endl << withoutZerosMatrix << endl;
			cout << "Matrix de canal" << endl << withZerosMatrix << endl;
		#endif

		nextNotUsedColumn += antennasChannelOrders[iAntenna];

		#ifdef DEBUG
			cout << "next....: " << nextNotUsedColumn << endl;
			getchar();
		#endif
	}
}
