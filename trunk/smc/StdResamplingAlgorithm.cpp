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
#include "StdResamplingAlgorithm.h"

using namespace std;

void StdResamplingAlgorithm::Resampling(ParticleWithChannelEstimation ***particles,int nParticles,std::vector<int> resamplingIndexes,std::vector<int> indexes)
{
	if(resamplingIndexes.size()!=indexes.size())
		throw RuntimeException("StdResamplingAlgorithm::Resampling: the size of the indexes vector and resampling indexes vector don't match.");

	int nParticlesToBeResampled = indexes.size();

	ParticleWithChannelEstimation **resParticles = new ParticleWithChannelEstimation*[nParticles];

	// the particles given by indexes are resampled
	for(int iParticle=0;iParticle<nParticlesToBeResampled;iParticle++)
	{
		resParticles[indexes[iParticle]] = ((*particles)[resamplingIndexes[iParticle]])->Clone();
		resParticles[indexes[iParticle]]->SetWeight(1.0/(double)nParticlesToBeResampled);
	}

	// the particles out of index are left the same. Their memory will not be released later
	int previousResampledParticle = 0;
	for(int iParticle=0;iParticle<nParticlesToBeResampled;iParticle++)
	{
		while(previousResampledParticle<indexes[iParticle])
		{
			resParticles[previousResampledParticle] = (*particles)[previousResampledParticle];
			previousResampledParticle++;
		}
		previousResampledParticle++;
	}
	while(previousResampledParticle<nParticles)
	{
		resParticles[previousResampledParticle] = (*particles)[previousResampledParticle];
		previousResampledParticle++;		
	}

	// the memory of the particles given by index is released
	for(int iParticle=0;iParticle<nParticlesToBeResampled;iParticle++)
		delete (*particles)[indexes[iParticle]];

	delete[] (*particles);
	*particles = resParticles;
}

void StdResamplingAlgorithm::Resampling(ParticleWithChannelEstimation ***particles,int nParticles,vector<int> indexes)
{
        ParticleWithChannelEstimation **resParticles = new ParticleWithChannelEstimation*[nParticles];

        for(int iParticle=0;iParticle<nParticles;iParticle++)
        {
                resParticles[iParticle] = ((*particles)[indexes[iParticle]])->Clone();
                resParticles[iParticle]->SetWeight(1.0/(double)nParticles);
        }

        for(int iParticle=0;iParticle<nParticles;iParticle++)
                delete (*particles)[iParticle];

        delete[] (*particles);
        *particles = resParticles;
}

