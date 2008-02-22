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

ParticleFilter::ParticleFilter(int nParticles):_capacity(nParticles),_nParticles(0),_particles(new ParticleWithChannelEstimation*[nParticles])
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

void ParticleFilter::Clear()
{
    for(uint i=0;i<_nParticles;i++)
    {
        delete _particles[i];
        _particles[i] = NULL;
    }
	_nParticles = 0;
}

void ParticleFilter::KeepParticles(std::vector<int> resamplingIndexes,std::vector<int> indexes)
{
	if(resamplingIndexes.size()!=indexes.size())
		throw RuntimeException("ParticleFilter::KeepParticles: the size of the indexes vector and resampling indexes vector don't match.");

	int nParticlesToBeResampled = indexes.size();

	ParticleWithChannelEstimation **resParticles = new ParticleWithChannelEstimation*[_capacity];

	// the particles given by indexes are resampled
	for(int iParticle=0;iParticle<nParticlesToBeResampled;iParticle++)
	{
		resParticles[indexes[iParticle]] = (_particles[resamplingIndexes[iParticle]])->Clone();
		resParticles[indexes[iParticle]]->SetWeight(1.0/(double)nParticlesToBeResampled);
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

void ParticleFilter::KeepParticles(vector<int> indexes)
{
	if(indexes.size()>_capacity)
		throw RuntimeException("ParticleFilter::KeepParticles: the number of selected particles is bigger than the number of particles in the filter.");

	#ifdef DEBUG
		cout << "indexes.size() = " << indexes.size() << " _nParticles = " << _nParticles << " _capacity = " << _capacity << endl;
	#endif

	ParticleWithChannelEstimation **resParticles = new ParticleWithChannelEstimation*[_capacity];

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

	for(uint iParticle=0;iParticle<_nParticles;iParticle++)
			delete _particles[iParticle];

	delete[] _particles;
	_particles = resParticles;
	_nParticles = indexes.size();
}

int ParticleFilter::iBestParticle()
{
// 	int iBestParticle;
// 	Util::Max(GetWeightsVector(),iBestParticle);
// 	return iBestParticle;

#ifdef DEBUG3
	cout << "_nParticles = " << _nParticles << endl;
#endif

	vector<bool> particleAlreadyCounted(_nParticles,false);
	vector<double> accumulatedWeights(_nParticles,0.0);
	uint iParticle,iTestedParticle;
	for(iParticle=0;iParticle<_nParticles;iParticle++)
	{
#ifdef DEBUG3
// 		cout << "iParticle = " << iParticle << endl;
#endif
		if(particleAlreadyCounted[iParticle])
			continue;

		accumulatedWeights[iParticle] = _particles[iParticle]->GetWeight();

		for(iTestedParticle=iParticle+1;iTestedParticle<_nParticles;iTestedParticle++)
		{
			if(particleAlreadyCounted[iTestedParticle])
				continue;
#ifdef DEBUG3
			cout << GetParticle(iParticle)->GetAllSymbolVectors();
			cout << "----------------" << endl;
			cout << GetParticle(iTestedParticle)->GetAllSymbolVectors();
// 			cout << "Una tecla..."; getchar();
#endif
			if(GetParticle(iParticle)->GetAllSymbolVectors().equal_to(GetParticle(iTestedParticle)->GetAllSymbolVectors()))
			{
#ifdef DEBUG2
				cout << "La partícula " << iParticle << " y la " << iTestedParticle << " son iguales" << endl;
#endif
				accumulatedWeights[iParticle] += GetParticle(iTestedParticle)->GetWeight();
				particleAlreadyCounted[iTestedParticle] = true;
			}
		}
	}
	int iBestParticle = Util::Max(accumulatedWeights);
#ifdef DEBUG2
	cout << "iBestParticle " << iBestParticle << " iBestParticle2 " << iBestParticle2 << endl;
	cout << "Los pesos: " << endl << GetWeightsVector();
	cout << "Los pesos acumulados" << endl;
	Util::Print(accumulatedWeights);
	getchar();
#endif
	return iBestParticle;
}
