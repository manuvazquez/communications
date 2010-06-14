/*
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License version 2 as published by the Free Software Foundation.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public License
   along with this library; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/

#ifndef PARTICLEWITHMULTIPLECHANNELSESTIMATIONANDMULTIPLECHANNELORDERAPP_H
#define PARTICLEWITHMULTIPLECHANNELSESTIMATIONANDMULTIPLECHANNELORDERAPP_H

#include <Particle.h>
#include <WithMultipleChannelsEstimationParticleAddon.h>
#include <WithMultipleChannelOrderAppParticleAddon.h>


class ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp : public Particle, public WithMultipleChannelsEstimationParticleAddon, public WithMultipleChannelOrderAppParticleAddon
{
public:
  ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp(double weight, int symbolVectorLength, int nTimeInstants, std::vector<std::vector <ChannelMatrixEstimator *> > channelMatrixEstimators);
  
  ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp* clone()
  {
	return new ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp(*this);
  }
};

#endif // PARTICLEWITHMULTIPLECHANNELSESTIMATIONANDMULTIPLECHANNELORDERAPP_H
