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
#include "WithChannelEstimationParticleAddon.h"

WithChannelEstimationParticleAddon::WithChannelEstimationParticleAddon(ChannelMatrixEstimator *channelMatrixEstimator,uint trajectorylength):_channelMatrixEstimators(1),_estimatedChannelMatrices(1,vector<MatrixXd>(trajectorylength))
#ifdef SAVE_CHANNEL_ESTIMATES_VARIANCES
,
_channelEstimatesVariances(1,vector<MatrixXd>(trajectorylength))
#endif
{
    _channelMatrixEstimators[0] = channelMatrixEstimator;
}

WithChannelEstimationParticleAddon::WithChannelEstimationParticleAddon(std::vector <ChannelMatrixEstimator *> channelMatrixEstimators,uint trajectorylength):_channelMatrixEstimators(channelMatrixEstimators),_estimatedChannelMatrices(channelMatrixEstimators.size(),vector<MatrixXd>(trajectorylength))
#ifdef SAVE_CHANNEL_ESTIMATES_VARIANCES
,
_channelEstimatesVariances(channelMatrixEstimators.size(),vector<MatrixXd>(trajectorylength))
#endif
{
}

WithChannelEstimationParticleAddon::WithChannelEstimationParticleAddon(const WithChannelEstimationParticleAddon& withChannelEstimationParticleAddon):_channelMatrixEstimators(withChannelEstimationParticleAddon._channelMatrixEstimators.size()),_estimatedChannelMatrices(withChannelEstimationParticleAddon._estimatedChannelMatrices)
#ifdef SAVE_CHANNEL_ESTIMATES_VARIANCES
,
_channelEstimatesVariances(withChannelEstimationParticleAddon._channelEstimatesVariances)
#endif
{
    for(uint i=0;i<_channelMatrixEstimators.size();i++)
        _channelMatrixEstimators[i] = withChannelEstimationParticleAddon._channelMatrixEstimators[i]->clone();
}

WithChannelEstimationParticleAddon::~WithChannelEstimationParticleAddon()
{
    for(uint i=0;i<_channelMatrixEstimators.size();i++)
        delete _channelMatrixEstimators[i];   
}
