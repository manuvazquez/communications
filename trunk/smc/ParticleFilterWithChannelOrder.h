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
#ifndef PARTICLEFILTERWITHCHANNELORDER_H
#define PARTICLEFILTERWITHCHANNELORDER_H

#include <ParticleFilter.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <ParticleWithChannelEstimationAndChannelOrder.h>

class ParticleFilterWithChannelOrder : public ParticleFilter
{
protected:
    vector<int> _candidateOrders;
    int _maxOrder,*_channelOrder2index,*_nParticlesPerChannelOrder;
    double *_channelOrderWeightsSum;
public:
    ParticleFilterWithChannelOrder(int nParticles,vector<int> candidateOrders);
    ~ParticleFilterWithChannelOrder();

    void SetParticle(ParticleWithChannelEstimation *particle,int n)
    {
        ParticleWithChannelEstimationAndChannelOrder *processedParticle;

        if(_particles[n]!=NULL)
        {
            // the counter for the channel order of the old particle is decreased
            processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particles[n]);
            _nParticlesPerChannelOrder[_channelOrder2index[processedParticle->GetChannelOrder()]]--;

            delete _particles[n];
        }

        // the counter for the channel order of the new particle is increased
        processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (particle);
        _nParticlesPerChannelOrder[_channelOrder2index[processedParticle->GetChannelOrder()]]++;

        _particles[n] = particle;

    }

    int NparticlesOfChannelOrderIndex(int iChannelOrder) { return _nParticlesPerChannelOrder[iChannelOrder];}

    int NchannelOrders() {return _candidateOrders.size();}

    vector<vector<int> > GetIndexesOfChannelOrders();
    void KeepParticles(std::vector<int> resamplingIndexes,std::vector<int> indexes);
    void KeepParticles(std::vector<int> resamplingIndexes);
//     void NormalizeParticlesByOrder();

};

#endif