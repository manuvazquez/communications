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

    _nParticlesPerChannelOrder = new int[_candidateOrders.size()];
    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _nParticlesPerChannelOrder[iChannelOrder] = 0;

    _channelOrderWeightsSum = new double[_maxOrder+1];
}


ParticleFilterWithChannelOrder::~ParticleFilterWithChannelOrder()
{
    delete[] _channelOrder2index;
    delete[] _nParticlesPerChannelOrder;
    delete[] _channelOrderWeightsSum;
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

void ParticleFilterWithChannelOrder::KeepParticles(std::vector<int> resamplingIndexes,std::vector<int> indexes)
{
    ParticleWithChannelEstimationAndChannelOrder *processedParticle;

    for(int iParticle=0;iParticle<resamplingIndexes.size();iParticle++)
    {
        // particle to be replaced
        processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particles[indexes[iParticle]]);

        _nParticlesPerChannelOrder[_channelOrder2index[processedParticle->GetChannelOrder()]]--;

        // replacing particle
        processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particles[resamplingIndexes[iParticle]]);
        _nParticlesPerChannelOrder[_channelOrder2index[processedParticle->GetChannelOrder()]]++;
    }

    ParticleFilter::KeepParticles(resamplingIndexes,indexes);
}

void ParticleFilterWithChannelOrder::KeepParticles(std::vector<int> resamplingIndexes)
{
    ParticleWithChannelEstimationAndChannelOrder *processedParticle;

    for(int i=0;i<_candidateOrders.size();i++)
        _nParticlesPerChannelOrder[_channelOrder2index[_candidateOrders[i]]] = 0;

    for(int iParticle=0;iParticle<_nParticles;iParticle++)
    {
        processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particles[resamplingIndexes[iParticle]]);
        _nParticlesPerChannelOrder[_channelOrder2index[processedParticle->GetChannelOrder()]]++;
    }

    ParticleFilter::KeepParticles(resamplingIndexes);
}

// void ParticleFilterWithChannelOrder::NormalizeParticlesByOrder()
// {
//     int iParticle;
//
//     // the sum of the weights of each group of particles is set to zero...
//     for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
//         _channelOrderWeightsSum[_candidateOrders[iChannelOrder]] = 0.0;
//
//     // ... to compute it here
//     for(iParticle=0;iParticle<_nParticles;iParticle++)
//     {
//         ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particles[iParticle]);
//
//         _channelOrderWeightsSum[processedParticle->GetChannelOrder()] += processedParticle->GetWeight();
//     }
//
//     // each particle is normalized according to its channel order
//     for(iParticle=0;iParticle<_nParticles;iParticle++)
//     {
//         ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particles[iParticle]);
//
//         processedParticle->SetWeight(processedParticle->GetWeight()/_channelOrderWeightsSum[processedParticle->GetChannelOrder()]);
//     }
// }
