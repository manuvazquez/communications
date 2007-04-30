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
#include "PSPBasedSMCAlgorithm.h"

PSPBasedSMCAlgorithm::PSPBasedSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances): SMCAlgorithm(name, alphabet, L, N, K, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances)
{
}


PSPBasedSMCAlgorithm::~PSPBasedSMCAlgorithm()
{
}


void PSPBasedSMCAlgorithm::InitializeParticles()
{
	// we begin with only one particle
	_particleFilter->SetParticle(new ParticleWithChannelEstimation(1.0,_N,_K,_channelEstimator->Clone()),0);
	_particleFilter->GetParticle(0)->SetSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);
}

void PSPBasedSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	uint nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);

	struct{
		int fromParticle;
		tVector symbolVector;
		double weight;
	}tParticleCandidate;

// 	vector<tParticleCandidate>

	// at first, there is only one particle
	int nParticles = 1;
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K+_d;iObservationToBeProcessed++)
	{
	}
}

