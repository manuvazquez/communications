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
#ifndef PARTICLEWITHCHANNELESTIMATIONANDLINEARDETECTIONANDACTIVEUSERS_H
#define PARTICLEWITHCHANNELESTIMATIONANDLINEARDETECTIONANDACTIVEUSERS_H

#include <Particle.h>
#include <WithChannelEstimationParticleAddon.h>
#include <WithLinearDetectionParticleAddon.h>
#include <WithActiveUsersParticleAddon.h>

/**
It implements a...

	@author Manu <manu@rustneversleeps>
*/
class ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers : public Particle, public WithChannelEstimationParticleAddon, public WithLinearDetectionParticleAddon, public WithActiveUsersParticleAddon
{
public:
    ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers(double weight, uint symbolVectorLength, uint nTimeInstants, ChannelMatrixEstimator* channelMatrixEstimator, LinearDetector* linearDetector);

    ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers* clone()
    {
        return new ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers(*this);
    }
};

#endif
