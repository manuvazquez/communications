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

// #define DEBUG2

EstimatedMIMOChannel::EstimatedMIMOChannel(int nInputs, int nOutputs, int memory, int length, int preambleLength, const ChannelMatrixEstimator *channelMatrixEstimator, const tMatrix &symbols, const tMatrix &observations, const vector<double> &noiseVariances): StillMemoryMIMOChannel(nInputs, nOutputs, memory, length),_channelMatrices(new tMatrix[length])
{
	int i;

	ChannelMatrixEstimator *channelMatrixEstimatorClone = channelMatrixEstimator->Clone();


	for(i=0;i<preambleLength;i++)
		_channelMatrices[i] = 0.0;


	tRange rAll,rInvolvedSymbolVectors(preambleLength-memory+1,preambleLength);
	for(i=preambleLength;i<_length;i++)
	{
		_channelMatrices[i] = channelMatrixEstimatorClone->nextMatrix(observations.col(i),symbols(rAll,rInvolvedSymbolVectors),noiseVariances[i]);
		rInvolvedSymbolVectors = rInvolvedSymbolVectors + 1;
#ifdef DEBUG2
		_channelMatrices[i](tRange(),tRange(0,1)) = 0;
#endif
#ifdef DEBUG
		cout << _channelMatrices[i] << endl;
#endif
	}

	delete channelMatrixEstimatorClone;
}


EstimatedMIMOChannel::~EstimatedMIMOChannel()
{
	delete[] _channelMatrices;
}


