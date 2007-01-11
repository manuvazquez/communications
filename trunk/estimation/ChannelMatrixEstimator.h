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
#ifndef CHANNELMATRIXESTIMATOR_H
#define CHANNELMATRIXESTIMATOR_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <types.h>

class ChannelMatrixEstimator{
protected:
	int _L,_Nm;
	tMatrix _lastEstimatedChannelMatrix;

	ChannelMatrixEstimator();
public:
	// initialEstimation is basically what LastEstimatedChannelMatrix is going to return when NextMatrix hasn't yet been called
    ChannelMatrixEstimator(tMatrix initialEstimation);
	virtual ~ChannelMatrixEstimator() {};

	virtual tMatrix NextMatrix(const tVector &observations,const tMatrix &symbolsMatrix,double noiseVariance) = 0;
	virtual ChannelMatrixEstimator *Clone() = 0;
	int Cols() { return _Nm;}
	int Rows() { return _L;}
	tMatrix LastEstimatedChannelMatrix() { return _lastEstimatedChannelMatrix;}

};

#endif
