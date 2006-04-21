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

void StdResamplingAlgorithm::Resampling(tMatrix  ***estimatedChannelMatrices,tMatrix ***detectedSymbols,ChannelMatrixEstimator ***particlesChannelMatrixEstimators,vector<int> indexes,int nParticles,int startResamplingTime,int endResamplingTime,int nTimeInstants)
{
	int iParticle,j;

	tMatrix **resEstimatedChannelMatrices;
	tMatrix **resDetectedSymbols;
	ChannelMatrixEstimator **resParticlesChannelMatrixEstimators;

	// new memory is allocated
	resEstimatedChannelMatrices = new tMatrix*[nParticles];
	resDetectedSymbols = new tMatrix*[nParticles];
	resParticlesChannelMatrixEstimators = new ChannelMatrixEstimator*[nParticles];

	for(iParticle=0;iParticle<nParticles;iParticle++)
	{

		// memory is allocated for each trajectory within a particle
		resEstimatedChannelMatrices[iParticle] = new tMatrix[nTimeInstants];

		// channel matrices
		for(j=startResamplingTime;j<endResamplingTime;j++)
		{
			resEstimatedChannelMatrices[iParticle][j] = (*estimatedChannelMatrices)[indexes[iParticle]][j];
		}

		// symbol vectors
		resDetectedSymbols[iParticle] = new tMatrix(*((*detectedSymbols)[indexes[iParticle]]));

		// channel matrix estimators
		resParticlesChannelMatrixEstimators[iParticle] = ((*particlesChannelMatrixEstimators)[indexes[iParticle]])->Clone();
	}

	// old memory is freed
	for(iParticle=0;iParticle<nParticles;iParticle++)
	{
		delete[] (*estimatedChannelMatrices)[iParticle];
		delete (*particlesChannelMatrixEstimators)[iParticle];
		delete (*detectedSymbols)[iParticle];
	}
	delete[] (*estimatedChannelMatrices);
	delete[] (*particlesChannelMatrixEstimators);
	delete[] (*detectedSymbols);

	*estimatedChannelMatrices = resEstimatedChannelMatrices;
	*detectedSymbols = resDetectedSymbols;
	*particlesChannelMatrixEstimators = resParticlesChannelMatrixEstimators;
}

