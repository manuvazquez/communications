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
#include "EstimatedMIMOChannel.h"

EstimatedMIMOChannel::EstimatedMIMOChannel(int nTx, int nRx, int memory, int length, const ChannelMatrixEstimator *channelMatrixEstimator, const tMatrix &symbols, const tMatrix &observations, const vector<double> &noiseVariances): StillMemoryMIMOChannel(nTx, nRx, memory, length),_channelMatrices(new tMatrix[length])
{
	int i;

	ChannelMatrixEstimator *channelMatrixEstimatorClone = channelMatrixEstimator->Clone();


	for(i=0;i<_memory-1;i++)
		_channelMatrices[i] = 0.0;


	tRange rAll,rInvolvedSymbolVectors(0,_memory-1);
	for(i=_memory-1;i<_length;i++)
	{
		_channelMatrices[i] = channelMatrixEstimatorClone->NextMatrix(observations.col(i),symbols(rAll,rInvolvedSymbolVectors),noiseVariances[i]);
		rInvolvedSymbolVectors = rInvolvedSymbolVectors + 1;
	}

	delete channelMatrixEstimatorClone;
}


EstimatedMIMOChannel::~EstimatedMIMOChannel()
{
	delete[] _channelMatrices;
}


