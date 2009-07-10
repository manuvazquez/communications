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
#include "MLSDmFeedBackAlgorithm.h"

MLSDmFeedBackAlgorithm::MLSDmFeedBackAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, double ARcoefficient, double samplingVariance, double ARprocessVariance): MLSDmAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm, ARcoefficient, samplingVariance, ARprocessVariance)
{
}

void MLSDmFeedBackAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
    MLSDmAlgorithm::Process(observations, noiseVariances);

	ParticleWithChannelEstimationAndChannelOrderAPP *bestParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->GetBestParticle()->clone());

	int nParticles = _particleFilter->Capacity();

	delete _particleFilter;

	for(int iObservationToBeProcessed=_iLastSymbolVectorToBeDetected+_d-2;iObservationToBeProcessed>=_startDetectionTime;iObservationToBeProcessed--)
	{
		for(int iChannelOrder=0;iChannelOrder<bestParticle->nChannelMatrixEstimators();iChannelOrder++)
		{
			tMatrix symbolVectors = bestParticle->GetSymbolVectors(tRange(iObservationToBeProcessed-_candidateOrders[iChannelOrder]+1,iObservationToBeProcessed));
			bestParticle->getChannelMatrixEstimator(iChannelOrder)->nextMatrix(observations.col(iObservationToBeProcessed),symbolVectors,noiseVariances[iObservationToBeProcessed]);
		}
	}

	bestParticle->SetWeight(1.0);

// 	// the available APP's just before the _startDetectionTime instant are copied into the particle
// 	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
// 		bestParticle->setChannelOrderAPP(_channelOrderAPPs(iChannelOrder,_startDetectionTime-1),iChannelOrder);

	_particleFilter = new ParticleFilter(nParticles);

	_particleFilter->AddParticle(bestParticle);

    MLSDmAlgorithm::Process(observations, noiseVariances);
}

