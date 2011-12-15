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
#ifndef ESTIMATEDMIMOCHANNEL_H
#define ESTIMATEDMIMOCHANNEL_H

#include <StillMemoryMIMOChannel.h>
#include <ChannelMatrixEstimator.h>

/**
    @author Manu <manu@rustneversleeps>
*/
class EstimatedMIMOChannel : public StillMemoryMIMOChannel
{
protected:
    std::vector<MatrixXd> _channelMatrices;
public:
    EstimatedMIMOChannel(uint nInputs, uint nOutputs, uint memory, uint length, uint preambleLength, const ChannelMatrixEstimator *channelMatrixEstimator, const MatrixXd &symbols, const MatrixXd &observations, const vector<double> &noiseVariances);
	
	virtual std::string name() const { return string("Estimated channel"); }

    MatrixXd at(uint n) const { return _channelMatrices[n];}
};

#endif
