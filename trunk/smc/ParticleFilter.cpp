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

void ParticleFilter::SelectParticles(std::vector<int> resamplingIndexes,std::vector<int> indexes)
{
	if(resamplingIndexes.size()!=indexes.size())
		throw RuntimeException("StdResamplingAlgorithm::Resampling: the size of the indexes vector and resampling indexes vector don't match.");

	int nParticlesToBeResampled = indexes.size();

	ParticleWithChannelEstimation **resParticles = new ParticleWithChannelEstimation*[_nParticles];

	// the particles given by indexes are resampled
	for(int iParticle=0;iParticle<nParticlesToBeResampled;iParticle++)
	{
		resParticles[indexes[iParticle]] = (_particles[resamplingIndexes[iParticle]])->Clone();
		resParticles[indexes[iParticle]]->SetWeight(1.0/(double)nParticlesToBeResampled);
	}

	// the particles out of index are left the same. Their memory will not be released later
	int previousResampledParticle = 0;
	for(int iParticle=0;iParticle<nParticlesToBeResampled;iParticle++)
	{
		while(previousResampledParticle<indexes[iParticle])
		{
			resParticles[previousResampledParticle] = _particles[previousResampledParticle];
			previousResampledParticle++;
		}
		previousResampledParticle++;
	}
	while(previousResampledParticle<_nParticles)
	{
		resParticles[previousResampledParticle] = _particles[previousResampledParticle];
		previousResampledParticle++;		
	}

	// the memory of the particles given by index is released
	for(int iParticle=0;iParticle<nParticlesToBeResampled;iParticle++)
		delete _particles[indexes[iParticle]];

	delete[] _particles;
	_particles = resParticles;
}

void ParticleFilter::SelectParticles(vector<int> indexes)
{
        ParticleWithChannelEstimation **resParticles = new ParticleWithChannelEstimation*[_nParticles];

        for(int iParticle=0;iParticle<_nParticles;iParticle++)
        {
                resParticles[iParticle] = (_particles[indexes[iParticle]])->Clone();
                resParticles[iParticle]->SetWeight(1.0/(double)_nParticles);
        }

        for(int iParticle=0;iParticle<_nParticles;iParticle++)
                delete _particles[iParticle];

        delete[] _particles;
        _particles = resParticles;
}

// void ParticleFilter::Resampling()
// {
// //     if(_resamplingCriterion.ResamplingNeeded(_particles,_nParticles))
//     if(_resamplingCriterion.ResamplingNeeded(GetWeightsVector()))
//     {
//         vector<int> indexes = StatUtil::Discrete_rnd(_nParticles,GetWeightsVector());
// 		SelectParticles(indexes);
//     }
// }
