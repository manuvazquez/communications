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
#include "ARchannel.h"

// #define DEBUG2

using namespace std;

ARchannel::ARchannel(int nTx, int nRx, int memory, int length,tMatrix initializationMatrix,vector<double> ARcoefficients,double ARvariance): StillMemoryMIMOChannel(nTx, nRx, memory, length),
//ARprocess constructor call
_ARproc(initializationMatrix,ARcoefficients,ARvariance)
{
	if(initializationMatrix.rows()!=nRx || initializationMatrix.cols()!=(nTx*memory))
		throw RuntimeException("ARchannel::ARchannel: the initialization matrix dimensions are not coherent with the received parameters.");

	_channelMatrices = new tMatrix[length];

	//initialization
	for(int i=_memory-1;i<_length;i++)
			_channelMatrices[i] = _ARproc.NextMatrix();
}

ARchannel::ARchannel(int nTx, int nRx, int memory, int length,ARprocess ARproc): StillMemoryMIMOChannel(nTx, nRx, memory, length),_ARproc(ARproc)
{
#ifdef DEBUG
	cout << "Principio de ARchannel" << endl;
#endif
	if(ARproc.Rows()!=nRx || ARproc.Cols()!=(nTx*memory))
		throw RuntimeException("ARchannel::ARchannel: the passed AR process is not compatible with the dimensions of the channel.");

	_channelMatrices = new tMatrix[length];

#ifdef DEBUG
	cout << "Antes del for en ARchannel" << endl;
#endif

	//initialization
	for(int i=_memory-1;i<_length;i++)
	{
			_channelMatrices[i] = _ARproc.NextMatrix();
#ifdef DEBUG2
		cout << _channelMatrices[i] << endl;
		cout << "Una tecla..."; getchar();
#endif
	}

#ifdef DEBUG
	cout << "Fin de ARchannel" << endl;
	cout << "Generado el canal" << endl;
#endif
}

ARchannel::ARchannel(const ARchannel &archannel):StillMemoryMIMOChannel(archannel),_channelMatrices(new tMatrix[_length]),_ARproc(archannel._ARproc)
{
	for(int i=_memory-1;i<_length;i++)
			_channelMatrices[i] = archannel._channelMatrices[i];
}

ARchannel::~ ARchannel()
{
	delete[] _channelMatrices;
}
