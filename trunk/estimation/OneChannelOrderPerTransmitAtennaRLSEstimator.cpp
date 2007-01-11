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
#include "OneChannelOrderPerTransmitAtennaRLSEstimator.h"

OneChannelOrderPerTransmitAtennaRLSEstimator::OneChannelOrderPerTransmitAtennaRLSEstimator(const tMatrix& initialEstimation, double forgettingFactor,const vector<int> &antennasChannelOrders): RLSEstimator(forgettingFactor),_antennasChannelOrders(antennasChannelOrders)
{
	_L = initialEstimation.rows();
	_Nm = initialEstimation.cols();
	int maxOrder = Util::Max(_antennasChannelOrders);
	int N = _antennasChannelOrders.size();
	if((_Nm % N)!=0)
		throw RuntimeException("OneChannelOrderPerTransmitAtennaRLSEstimator::OneChannelOrderPerTransmitAtennaRLSEstimator: the length of the antennas channel orders vector is not coherent with the initial estimation channel matrix dimensions");
}


OneChannelOrderPerTransmitAtennaRLSEstimator::~OneChannelOrderPerTransmitAtennaRLSEstimator()
{
}


tMatrix OneChannelOrderPerTransmitAtennaRLSEstimator::NextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
    return RLSEstimator::NextMatrix(observations, symbolsMatrix, noiseVariance);
}

