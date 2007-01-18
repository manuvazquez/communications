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

OneChannelOrderPerTransmitAtennaMIMOChannel::OneChannelOrderPerTransmitAtennaMIMOChannel(int nTx, int nRx, int length,const std::vector<int> &antennasChannelOrders): StillMemoryMIMOChannel(nTx, nRx,antennasChannelOrders[Util::Max(antennasChannelOrders)],length)
,_antennasChannelOrders(antennasChannelOrders),_nChannelMatrixNotNullColumns(Util::Sum(antennasChannelOrders))
{
// 	_maxChannelOrder = -1;
// 	for(int i=0;i<_nTx;i++)
// 	{
// 		if(_antennasChannelOrders[i] > _maxChannelOrder)
// 			_maxChannelOrder = _antennasChannelOrders[i];
//
// 		_nChannelMatrixNotNullColumns += _antennasChannelOrders[i];
// 	}
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
		throw RuntimeException("OneChannelOrderPerTransmitAtennaMIMOChannel::WithoutZerosMatrixToWithZerosMatrix: resultant matrix has not enough columns");

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

tMatrix OneChannelOrderPerTransmitAtennaMIMOChannel::WithZerosMatrixToWithoutZerosMatrix(const tMatrix &withZerosMatrix,int N,const vector<int> &antennasChannelOrders)
{
        int maxChannelOrder = antennasChannelOrders[Util::Max(antennasChannelOrders)];
        int channelOrdersSum = Util::Sum(antennasChannelOrders);

        tMatrix withoutZerosMatrix(withZerosMatrix.rows(),channelOrdersSum);

        #ifdef DEBUG
                cout << "las columnas de la matriz origen: " << withZerosMatrix.cols() << endl;
                cout << "las columnas de la matriz destino: " << withoutZerosMatrix.cols() << endl;
                cout << "maxChannelOrder: " << maxChannelOrder << " N: " << N << endl;
        #endif

        tRange rAllRows(0,withoutZerosMatrix.rows()-1);

        withoutZerosMatrix = 0.0;

        int nextNotUsedColumn = 0;
        for(int iAntenna=0;iAntenna<N;iAntenna++)
        {
                tRange rSourceColumns((maxChannelOrder-antennasChannelOrders[iAntenna])*N+iAntenna,N*maxChannelOrder-1,N);

                #ifdef DEBUG
                        cout << "El rango " << endl << rSourceColumns << endl;
                        cout << withoutZerosMatrix(rAllRows,tRange(nextNotUsedColumn,nextNotUsedColumn+antennasChannelOrders[iAntenna]-1)) << endl;
                        cout << "hola" << endl;
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

        return withoutZerosMatrix;
}

void OneChannelOrderPerTransmitAtennaMIMOChannel::CompleteSymbolsVectorToOnlyInvolvedSymbolsVector(const tVector &withZerosSymbolsVector,int N,const vector<int> &antennasChannelOrders,tVector &withoutZerosSymbolsVector)
{
	int maxChannelOrder = antennasChannelOrders[Util::Max(antennasChannelOrders)];
	int channelOrdersSum = Util::Sum(antennasChannelOrders);

	if(withoutZerosSymbolsVector.size()!= channelOrdersSum)
		throw RuntimeException("ARoneChannelOrderPerTransmitAtennaMIMOChannel::CompleteSymbolsVectorToOnlyInvolvedSymbolsVector: resultant vector has not the right number of columns");

	withoutZerosSymbolsVector = 0.0;

	int nextNotUsedColumn = 0;
	for(int iAntenna=0;iAntenna<N;iAntenna++)
	{
		tRange rSourceColumns((maxChannelOrder-antennasChannelOrders[iAntenna])*N+iAntenna,N*maxChannelOrder-1,N);
		tRange rDestinyColumns(nextNotUsedColumn,nextNotUsedColumn+antennasChannelOrders[iAntenna]-1);

		withoutZerosSymbolsVector(rDestinyColumns).inject(withZerosSymbolsVector(rSourceColumns));

		nextNotUsedColumn += antennasChannelOrders[iAntenna];
	}
}
