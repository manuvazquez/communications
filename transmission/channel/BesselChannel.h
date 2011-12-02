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
#ifndef BESSELCHANNEL_H
#define BESSELCHANNEL_H

#include <StillMemoryMIMOChannel.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <StatUtil.h>
#include <DelayPowerProfile.h>
#include <ContinuousPowerProfile.h>
#include <Eigen/Cholesky>

class BesselChannel : public StillMemoryMIMOChannel
{
protected:
    std::vector<MatrixXd> _channelMatrices;   
public:
    BesselChannel(uint nInputs, uint nOutputs, uint memory, uint length, double velocity, double carrierFrequency, double T, const DelayPowerProfile &powerProfile);

    MatrixXd at(uint n) const { return _channelMatrices[n];}
    
    virtual void set(uint n, MatrixXd mat)
    {
	  if(mat.rows()!=_nOutputs || mat.cols()!=_nInputsMemory)
		throw RuntimeException("BesselChannel:set: matrix dimensions are wrong.");

	  _channelMatrices[n] = mat;
	}    
};

#endif
