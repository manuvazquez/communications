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
#include "ChannelMatrixEstimator.h"

ChannelMatrixEstimator::ChannelMatrixEstimator(int N):_L(0),_Nm(0),_N(N)
{
}

ChannelMatrixEstimator::ChannelMatrixEstimator(tMatrix initialEstimation,int N):_lastEstimatedChannelMatrix(initialEstimation),_L(initialEstimation.rows()),_Nm(initialEstimation.cols()),_N(N)
{
	if(_Nm < _N)
		throw RuntimeException("ChannelMatrixEstimator::ChannelMatrixEstimator: number of columns of \"initialEstimation\"  is less than N");
// 	if((_Nm % _N) != 0)
// 		throw RuntimeException("ChannelMatrixEstimator::ChannelMatrixEstimator: number of columns of \"initialEstimation\"  is not a multiple of N");
}


