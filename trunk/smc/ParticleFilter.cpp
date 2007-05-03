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

// #define DEBUG

ParticleFilter::ParticleFilter(int nParticles):_nParticles(nParticles),_nActualParticles(0),_particles(new ParticleWithChannelEstimation*[nParticles])
{
    for(int i=0;i<_nParticles;i++)
    {
        _particles[i] = NULL;
    }
}

ParticleFilter::~ParticleFilter()
{
    for(int i=0;i<_nActualParticles;i++)
    {
        delete _particles[i];
    }

    delete[] _particles;
}

void ParticleFilter::KeepParticles(std::vector<int> resamplingIndexes,std::vector<int> indexes)
{
	if(resamplingIndexes.size()!=indexes.size())
		throw RuntimeException("ParticleFilter::KeepParticles: the size of the indexes vector and resampling indexes vector don't match.");

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

void ParticleFilter::KeepParticles(vector<int> indexes)
{
	if(indexes.size()>_nParticles)
		throw RuntimeException("ParticleFilter::KeepParticles: the number of selected particles is bigger than the number of particles in the filter.");

	#ifdef DEBUG
		cout << "indexes.size() = " << indexes.size() << " _nActualParticles = " << _nActualParticles << " _nParticles = " << _nParticles << endl;
	#endif

	ParticleWithChannelEstimation **resParticles = new ParticleWithChannelEstimation*[_nParticles];

	vector<bool> particleNeeded(_nActualParticles,false);
	for(uint iParticle=0;iParticle<indexes.size();iParticle++)
		particleNeeded[indexes[iParticle]] = true;

	// the memory occupied by the particles that are not gonna be resampled is released
	for(int iParticle=0;iParticle<_nActualParticles;iParticle++)
		if(!particleNeeded[iParticle])
		{
			delete _particles[iParticle];
			_particles[iParticle] = NULL;
		}

	#ifdef DEBUG
		cout << "antes de replicar" << endl;
		cout << "El vector de bools es " << endl;
		Util::Print(particleNeeded);
	#endif

	for(uint iParticle=0;iParticle<indexes.size();iParticle++)
	{
			#ifdef DEBUG
				cout << "iParticle = " << iParticle << ". Accediendo a la particula " <<  indexes[iParticle] << endl;
			#endif
			resParticles[iParticle] = (_particles[indexes[iParticle]])->Clone();
			resParticles[iParticle]->SetWeight(1.0/(double)_nParticles);
	}

	#ifdef DEBUG
		cout << "despues de replicar" << endl;
	#endif

	for(int iParticle=0;iParticle<_nActualParticles;iParticle++)
			delete _particles[iParticle];

	delete[] _particles;
	_particles = resParticles;
	_nActualParticles = indexes.size();
}

// void ParticleFilter::KeepParticles(vector<int> indexes)
// {
// 	ParticleWithChannelEstimation **resParticles = new ParticleWithChannelEstimation*[_nParticles];
//
// 	vector<bool> particleNeeded(_nParticles,false);
// 	for(int iParticle=0;iParticle<_nParticles;iParticle++)
// 		particleNeeded[indexes[iParticle]] = true;
//
// 	// the memory occupied by the particles that are not gonna be resampled is released
// 	for(int iParticle=0;iParticle<_nParticles;iParticle++)
// 		if(!particleNeeded[iParticle])
// 		{
// 			delete _particles[iParticle];
// 			_particles[iParticle] = NULL;
// 		}
//
// 	for(int iParticle=0;iParticle<_nParticles;iParticle++)
// 	{
// 			resParticles[iParticle] = (_particles[indexes[iParticle]])->Clone();
// 			resParticles[iParticle]->SetWeight(1.0/(double)_nParticles);
// 	}
//
// 	for(int iParticle=0;iParticle<_nParticles;iParticle++)
// 			delete _particles[iParticle];
//
// 	delete[] _particles;
// 	_particles = resParticles;
// }
