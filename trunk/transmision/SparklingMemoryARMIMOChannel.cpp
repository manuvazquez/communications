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
#include "SparklingMemoryARMIMOChannel.h"

// the seed used to create the random objects is generated from the system time
// #define RANDOM_SEED

SparklingMemoryARMIMOChannel::SparklingMemoryARMIMOChannel(int nTx, int nRx, int length, std::vector< int > candidateOrders, int initialChannelOrderIndex,double mean,double variance,vector<double> ARcoefficients,double ARvariance): SparklingMemoryMIMOChannel(nTx, nRx, length, candidateOrders, initialChannelOrderIndex),_ARprocess(StatUtil::RandnMatrix(nRx,nTx*_maxOrder,mean,variance),
ARcoefficients,ARvariance)
{
	int presentChannelOrder,iPresentChannelOrder;
	tRange rAllObservationsRows(0,_nRx-1);
	double overAllNorm,lastTapNorm,ratio,sample;
	double changeProbability = 0.1;

	#ifdef RANDOM_SEED
		static Random randomGenerator;
	#else
		static Random randomGenerator(3);
	#endif

	_channelMatrices = new tMatrix[length];
	_channelOrders = new int[length];

	// the first channel matrix has memory "initialChannelOrderIndex"
	iPresentChannelOrder = initialChannelOrderIndex;
	presentChannelOrder = _candidateOrders[_initialChannelOrderIndex];

	for(int i=_maxOrder-1;i<_length;i++)
	{
		_channelOrders[i] = iPresentChannelOrder;

		cout << _candidateOrders[iPresentChannelOrder] << endl;

		_channelMatrices[i] = _ARprocess.NextMatrix()(rAllObservationsRows,tRange(0,_nTx*_candidateOrders[iPresentChannelOrder]-1));

		cout << "La matriz de canal" << endl << _channelMatrices[i] << endl;

		// overall Frobenius norm of the channel is computed
		overAllNorm = pow(Blas_NormF(_channelMatrices[i]),2.0);

		// Frobenius norm of the last "tap" of the channel
		lastTapNorm = pow(Blas_NormF(_channelMatrices[i](rAllObservationsRows,tRange((_candidateOrders[iPresentChannelOrder]-1)*nTx,_candidateOrders[iPresentChannelOrder]*nTx-1))),2.0);

		cout << overAllNorm << "," << lastTapNorm << endl;

		ratio = lastTapNorm/(overAllNorm/_candidateOrders[iPresentChannelOrder]);

		cout << "el ratio es " << ratio << endl;

		// if the energy of the last tap is relatively small
		if(ratio < 1.0)
		{
			if(_candidateOrders[iPresentChannelOrder] > 1)
				iPresentChannelOrder -= ((1.0 - ratio) > randomGenerator.rand()) && (randomGenerator.rand() < changeProbability);
		}
		// if the energy of the last tap is relatively big
		else
		{
			if(_candidateOrders[iPresentChannelOrder] < _maxOrder)
			{
				iPresentChannelOrder += ((ratio - 1.0) > randomGenerator.rand()) && (randomGenerator.rand() < changeProbability);
			}
		}

		getchar();
	}

	// the first channel order is modified in the loop but it must be "initialChannelOrderIndex"
	_channelOrders[_maxOrder-1] = initialChannelOrderIndex;
}

SparklingMemoryARMIMOChannel::SparklingMemoryARMIMOChannel(const SparklingMemoryARMIMOChannel &channel):SparklingMemoryMIMOChannel(channel),_ARprocess(channel._ARprocess),_channelMatrices(new tMatrix[_length]),_channelOrders(new int[_length])
{
	for(int i=_maxOrder-1;i<_length;i++)
	{
		_channelOrders[i] = channel._channelOrders[i];
		_channelMatrices[i] = channel._channelMatrices[i];
	}
}

SparklingMemoryARMIMOChannel::~SparklingMemoryARMIMOChannel()
{
	delete[] _channelMatrices;
	delete[] _channelOrders;
}

