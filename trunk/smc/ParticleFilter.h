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
#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <Particle.h>
#include <ParticleWithChannelEstimation.h>
#include <StdResamplingAlgorithm.h>
#include <ResamplingCriterion.h>
#include <StatUtil.h>

class ParticleFilter{
protected:
    int _nParticles;
    ResamplingCriterion _resamplingCriterion;
    StdResamplingAlgorithm _resamplingAlgorithm;
    ParticleWithChannelEstimation **_particles;
public:
    ParticleFilter(int nParticles,const ResamplingCriterion &resamplingCriterion,const StdResamplingAlgorithm &resamplingAlgorithm);

    ~ParticleFilter();

	void Resampling();
	ParticleWithChannelEstimation *GetParticle(int n) { return _particles[n];}

	void SetParticle(ParticleWithChannelEstimation *particle,int n)
	{
		delete _particles[n];
		_particles[n] = particle;
	}
	
    tVector GetWeightsVector() 
    {
        tVector weights(_nParticles);
        for(int i=0;i<_nParticles;i++)
            weights(i) = _particles[i]->GetWeight();
        return weights;
    }

    void NormalizeWeights()
    {
        double sum = 0.0;
        int i;

        for(i=0;i<_nParticles;i++)
            sum += _particles[i]->GetWeight();

        for(i=0;i<_nParticles;i++)
            _particles[i]->SetWeight(_particles[i]->GetWeight()/sum);
    }

	int Nparticles() { return _nParticles;}
};

#endif
