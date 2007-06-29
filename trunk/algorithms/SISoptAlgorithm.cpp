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
#include "SISoptAlgorithm.h"

SISoptAlgorithm::SISoptAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances): SMCAlgorithm(name, alphabet, L, N, K, m, channelEstimator, preamble, 0, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances)
{
}

void SISoptAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	int k,iParticle,iSampledVector;
	vector<tSymbol> testedVector(_N),sampledVector(_N);
	tRange mMinus1FirstColumns(0,_m-2);

	// it selects all rows in the symbols Matrix
	tRange rAllSymbolRows(0,_N-1);

	// it includes all symbol vectors involved in the smoothing
	tMatrix involvedSymbolVectors(_N,_m);

	uint nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);

	// a likelihood is computed for every possible symbol vector
	tVector likelihoods(nSymbolVectors);

	tRange mPrecedentColumns(_startDetectionTime-_m+1,_startDetectionTime);
	tRange mMinus1PrecedentColumns(_startDetectionTime-_m+1,_startDetectionTime-1);

	// for each time instant
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
		for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)
		{
			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			// the m-1 already detected symbol vectors are copied into the matrix:
			involvedSymbolVectors(rAllSymbolRows,mMinus1FirstColumns).inject(processedParticle->GetSymbolVectors(mMinus1PrecedentColumns));

			for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.IntToSymbolsArray(iTestedVector,testedVector);

				// current tested vector is copied in the m-th position
				for(k=0;k<_N;k++)
					involvedSymbolVectors(k,_m-1) = testedVector[k];

				likelihoods(iTestedVector) = processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->Likelihood(observations.col(iObservationToBeProcessed),involvedSymbolVectors,noiseVariances[iObservationToBeProcessed]);
			} // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

			tVector probabilities(nSymbolVectors);
			try {
				// probabilities are computed by normalizing the likelihoods
				probabilities = Util::Normalize(likelihoods);
			} catch (AllElementsNullException) {
				// if all the likelihoods are null
				probabilities = 1.0/(double)nSymbolVectors;
			}

			// one sample from the discrete distribution is taken
			iSampledVector = StatUtil::Discrete_rnd(probabilities);

			// the above index is turned into a vector
			_alphabet.IntToSymbolsArray(iSampledVector,sampledVector);

			// sampled symbols are copied into the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledVector);

			// channel matrix is estimated by means of the particle channel estimator
			processedParticle->SetChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->NextMatrix(observations.col(iObservationToBeProcessed),processedParticle->GetSymbolVectors(mPrecedentColumns),noiseVariances[iObservationToBeProcessed]));

			processedParticle->SetWeight(processedParticle->GetWeight()* Util::Sum(likelihoods));
		} // for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)

		mPrecedentColumns = mPrecedentColumns + 1;
		mMinus1PrecedentColumns = mMinus1PrecedentColumns + 1;

		_particleFilter->NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
// 			Resampling();
            _resamplingAlgorithm->ResampleWhenNecessary(_particleFilter);

	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
}

