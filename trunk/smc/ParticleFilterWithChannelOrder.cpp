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
#include "ParticleFilterWithChannelOrder.h"

ParticleFilterWithChannelOrder::ParticleFilterWithChannelOrder(int nParticles,vector<int> candidateOrders)
 : ParticleFilter(nParticles),_candidateOrders(candidateOrders),_maxOrder(-1)
{
    // it finds out the maximum channel order
    for(int i=0;i<_candidateOrders.size();i++)
        if(_candidateOrders[i]>_maxOrder)
            _maxOrder = _candidateOrders[i];

    // a vector that associate a channel order with its corresponding index is generated
    _channelOrder2index = new int[_maxOrder+1];
    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _channelOrder2index[_candidateOrders[iChannelOrder]] = iChannelOrder;
}


ParticleFilterWithChannelOrder::~ParticleFilterWithChannelOrder()
{
    delete[] _channelOrder2index;
}

/**
 * It returns a vector of vectors of int, such that vector[i] contains the indexes of the particles of order _candidateOrders[i]
 * @return 
 */
vector<vector<int> > ParticleFilterWithChannelOrder::GetIndexesOfChannelOrders()
{
    vector<vector<int> > res(_candidateOrders.size());

    for(int iParticle=0;iParticle<_nParticles;iParticle++)
    {
        ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particles[iParticle]);

        res[_channelOrder2index[processedParticle->GetChannelOrder()]].push_back(iParticle);
    }
    return res;
}
