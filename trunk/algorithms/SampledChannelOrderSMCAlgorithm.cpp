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
#include "SampledChannelOrderSMCAlgorithm.h"

SampledChannelOrderSMCAlgorithm::SampledChannelOrderSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_particleFilter(nParticles,_candidateOrders)
{
}

void SampledChannelOrderSMCAlgorithm::InitializeParticles()
{
    int iParticlePresentOrder,nParticlesPresentOrder,iParticle=0;
    int iChannelMatrixEstimator;

    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
    {
    // If the number of particles is not divisible by the number of channel orders, the remaining particles are assigned to the first channel orders (the + (......) term)
        nParticlesPresentOrder = _particleFilter.Nparticles()/_candidateOrders.size() + (_particleFilter.Nparticles() % _candidateOrders.size() > iChannelOrder);

        for(iParticlePresentOrder=0;iParticlePresentOrder<nParticlesPresentOrder;iParticlePresentOrder++,iParticle++)
        {
            // a clone of each of the channel matrix estimators is constructed...
            vector< ChannelMatrixEstimator * > thisParticleChannelMatrixEstimators(_candidateOrders.size());
            for(iChannelMatrixEstimator=0;iChannelMatrixEstimator<_candidateOrders.size();iChannelMatrixEstimator++)
                thisParticleChannelMatrixEstimators[iChannelMatrixEstimator] = _channelEstimators[iChannelMatrixEstimator]->Clone();

            // ... and passed within a vector to each particle
            _particleFilter.SetParticle(new ParticleWithChannelEstimationAndChannelOrder(1.0/(double)_particleFilter.Nparticles(),_N,_K,thisParticleChannelMatrixEstimators,_candidateOrders[iChannelOrder]),iParticle);
        }
    }
}

void SampledChannelOrderSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
}

