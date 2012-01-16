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

EstimatedMIMOChannel::EstimatedMIMOChannel(uint nInputs, uint nOutputs, uint memory, uint length, uint preambleLength, const ChannelMatrixEstimator *channelMatrixEstimator, const MatrixXd &symbols, const MatrixXd &observations, const vector<double> &noiseVariances): StillMemoryMIMOChannel(nInputs, nOutputs, memory, length),_channelMatrices(length)
{
	uint i;

	ChannelMatrixEstimator *channelMatrixEstimatorClone = channelMatrixEstimator->clone();


    MatrixXd nullMatrix = MatrixXd::Zero(1,1);
    nullMatrix.resize(0,0);
    
	for(i=0;i<preambleLength;i++)
		_channelMatrices[i] = nullMatrix;

//     MatrixXd symbols = Util::lapack2eigen(symbols);
//     MatrixXd observations = Util::lapack2eigen(observations);

	for(i=preambleLength;i<_length;i++)
        _channelMatrices[i] = channelMatrixEstimatorClone->nextMatrix(observations.col(i),symbols.block(0,i-memory+1,nInputs,memory),noiseVariances[i]);      

	delete channelMatrixEstimatorClone;
}
