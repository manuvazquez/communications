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
#include "ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation.h"

ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation::ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation(double weight, int symbolVectorLength, int nTimeInstants, std::vector< ChannelMatrixEstimator * > channelMatrixEstimators, std::vector< LinearDetector * > linearDetectors, ChannelOrderEstimator* channelOrderEstimator): ParticleWithChannelEstimationAndLinearDetection(weight, symbolVectorLength, nTimeInstants, channelMatrixEstimators, linearDetectors), WithChannelOrderEstimationParticleAddon(channelOrderEstimator)
{
}

ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation::ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation(const ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation& particle):ParticleWithChannelEstimationAndLinearDetection(particle),WithChannelOrderEstimationParticleAddon(particle)
{
}

ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation* ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation::clone()
{
	return new ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation(*this);
}
