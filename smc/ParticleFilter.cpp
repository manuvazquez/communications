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

// #define DEBUG3

ParticleFilter::ParticleFilter(uint nParticles):_capacity(nParticles),_nParticles(0),_particles(new Particle*[nParticles])
{
    for(uint i=0;i<_capacity;i++)
    {
        _particles[i] = NULL;
    }
}

ParticleFilter::~ParticleFilter()
{
    for(uint i=0;i<_nParticles;i++)
    {
        delete _particles[i];
    }

    delete[] _particles;
}

void ParticleFilter::clear()
{
    for(uint i=0;i<_nParticles;i++)
    {
        delete _particles[i];
        _particles[i] = NULL;
    }
    _nParticles = 0;
}

void ParticleFilter::keepParticles(std::vector<uint> resamplingIndexes,std::vector<uint> indexes)
{
    if(resamplingIndexes.size()!=indexes.size())
        throw RuntimeException("ParticleFilter::KeepParticles: the size of the indexes vector and resampling indexes vector don't match.");

    int nParticlesToBeResampled = indexes.size();

    Particle **resParticles = new Particle*[_capacity];

    // the particles given by indexes are resampled
    for(int iParticle=0;iParticle<nParticlesToBeResampled;iParticle++)
    {
        resParticles[indexes[iParticle]] = (_particles[resamplingIndexes[iParticle]])->clone();
        resParticles[indexes[iParticle]]->setWeight(1.0/(double)nParticlesToBeResampled);
    }

    // the particles out of index are left the same. Their memory will not be released later
    uint previousResampledParticle = 0;
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

void ParticleFilter::keepParticles(vector<uint> indexes)
{
    if(indexes.size()>_capacity)
        throw RuntimeException("ParticleFilter::KeepParticles: the number of selected particles is bigger than the number of particles in the filter.");

    Particle **resParticles = new Particle*[_capacity];

    vector<bool> particleNeeded(_nParticles,false);
    for(uint iParticle=0;iParticle<indexes.size();iParticle++)
        particleNeeded[indexes[iParticle]] = true;

    // the memory occupied by the particles that are not gonna be resampled is released
    for(uint iParticle=0;iParticle<_nParticles;iParticle++)
        if(!particleNeeded[iParticle])
        {
            delete _particles[iParticle];
            _particles[iParticle] = NULL;
        }

    for(uint iParticle=0;iParticle<indexes.size();iParticle++)
    {
            resParticles[iParticle] = (_particles[indexes[iParticle]])->clone();
            resParticles[iParticle]->setWeight(1.0/(double)_nParticles);
    }

    for(uint iParticle=0;iParticle<_nParticles;iParticle++)
            delete _particles[iParticle];

    delete[] _particles;
    _particles = resParticles;
    _nParticles = indexes.size();
}

int ParticleFilter::iBestParticle() const
{
//  int iBestParticle;
//  Util::max(getWeightsVector(),iBestParticle);
//  return iBestParticle;

    vector<bool> particleAlreadyCounted(_nParticles,false);
    vector<double> accumulatedWeights(_nParticles,0.0);
    uint iParticle,iTestedParticle;
    for(iParticle=0;iParticle<_nParticles;iParticle++)
    {
        if(particleAlreadyCounted[iParticle])
            continue;

        accumulatedWeights[iParticle] = _particles[iParticle]->getWeight();

        for(iTestedParticle=iParticle+1;iTestedParticle<_nParticles;iTestedParticle++)
        {
            if(particleAlreadyCounted[iTestedParticle])
                continue;

//             if(getParticle(iParticle)->getAllSymbolVectors().equal_to(getParticle(iTestedParticle)->getAllSymbolVectors()))
            if(getParticle(iParticle)->getSymbolVectors() == getParticle(iTestedParticle)->getSymbolVectors())
            {
                accumulatedWeights[iParticle] += getParticle(iTestedParticle)->getWeight();
                particleAlreadyCounted[iTestedParticle] = true;
            }
        }
    }
    int iBestParticle = Util::max(accumulatedWeights);

    return iBestParticle;
}

void ParticleFilter::printWeights() const
{
    VectorXd weights = getWeightsVector();
    
    for(uint i=0;i<_nParticles;i++)
        cout << i << "\t" << weights(i) << endl;
        
}
