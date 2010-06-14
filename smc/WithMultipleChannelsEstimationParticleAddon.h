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

#ifndef WITHMULTIPLECHANNELSESTIMATIONPARTICLEADDON_H
#define WITHMULTIPLECHANNELSESTIMATIONPARTICLEADDON_H

#include <WithChannelEstimationParticleAddon.h>
#include <ChannelMatrixEstimator.h>
#include <vector>

class WithMultipleChannelsEstimationParticleAddon
{
protected:
  std::vector<WithChannelEstimationParticleAddon> _channelEstimationParticleAddons;
public:
  WithMultipleChannelsEstimationParticleAddon(std::vector<std::vector <ChannelMatrixEstimator *> > channelMatrixEstimators,uint trajectoryLength);

    MatrixXd getChannelMatrix(uint iChannel,int iChannelOrder,int n) const
    {
        return _channelEstimationParticleAddons[iChannel].getChannelMatrix(iChannelOrder,n);
    }

    void setChannelMatrix(uint iChannel,int iChannelOrder,int n,const MatrixXd &matrix)
    {
	  _channelEstimationParticleAddons[iChannel].setChannelMatrix(iChannelOrder,n,matrix);
    }

    ChannelMatrixEstimator *getChannelMatrixEstimator(uint iChannel,int iChannelOrder) const
    { 
		return _channelEstimationParticleAddons[iChannel].getChannelMatrixEstimator(iChannelOrder);
    }
};

#endif // WITHMULTIPLECHANNELSESTIMATIONPARTICLEADDON_H
