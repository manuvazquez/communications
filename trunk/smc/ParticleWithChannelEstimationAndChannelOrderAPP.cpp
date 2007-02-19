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
#include "ParticleWithChannelEstimationAndChannelOrderAPP.h"

using namespace std;

ParticleWithChannelEstimationAndChannelOrderAPP::ParticleWithChannelEstimationAndChannelOrderAPP(double weight, int symbolVectorLength, int nTimeInstants, std::vector< ChannelMatrixEstimator * > channelMatrixEstimators): ParticleWithChannelEstimation(weight, symbolVectorLength, nTimeInstants, channelMatrixEstimators),WithChannelOrderAppParticleAddon(channelMatrixEstimators.size())
{
	for(uint iChannelOrder=0;iChannelOrder<_channelMatrixEstimators.size();iChannelOrder++)
		_channelOrderAPP[iChannelOrder] = 1.0/(double)_channelMatrixEstimators.size();
}

ParticleWithChannelEstimationAndChannelOrderAPP *ParticleWithChannelEstimationAndChannelOrderAPP::Clone()
{
	return new ParticleWithChannelEstimationAndChannelOrderAPP(*this);
}

ParticleWithChannelEstimationAndChannelOrderAPP::ParticleWithChannelEstimationAndChannelOrderAPP(const ParticleWithChannelEstimationAndChannelOrderAPP& particle):ParticleWithChannelEstimation(particle),WithChannelOrderAppParticleAddon(particle)
{
}
