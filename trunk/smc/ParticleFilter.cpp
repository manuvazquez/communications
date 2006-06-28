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
#include "ParticleFilter.h"

ParticleFilter::ParticleFilter(int nParticles,const ResamplingCriterion &resamplingCriterion,const StdResamplingAlgorithm &resamplingAlgorithm):_nParticles(nParticles),_resamplingCriterion(resamplingCriterion),_resamplingAlgorithm(resamplingAlgorithm),_particles(new ParticleWithChannelEstimation*[nParticles])
{
    for(int i=0;i<_nParticles;i++)
    {
        _particles[i] = NULL;
    }
}

ParticleFilter::~ParticleFilter()
{
    for(int i=0;i<_nParticles;i++)
    {
        delete _particles[i];
    }

    delete[] _particles;
}

void ParticleFilter::Resampling()
{
    if(_resamplingCriterion.ResamplingNeeded(_particles,_nParticles))
    {
        vector<int> indexes = StatUtil::Discrete_rnd(_nParticles,GetWeightsVector());

        _resamplingAlgorithm.Resampling(&_particles,_nParticles,indexes);
    }
}
