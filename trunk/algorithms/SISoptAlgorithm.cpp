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

SISoptAlgorithm::SISoptAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, 0, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances)
{
}

void SISoptAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	int k,iParticle,iSampledVector;
	vector<tSymbol> testedVector(_nInputs),sampledVector(_nInputs);
	tRange mMinus1FirstColumns(0,_channelOrder-2);

	// it selects all rows in the symbols Matrix
	tRange rAllSymbolRows(0,_nInputs-1);

	// it includes all symbol vectors involved in the smoothing
	tMatrix involvedSymbolVectors(_nInputs,_channelOrder);

	uint nSymbolVectors = (int) pow((double)_alphabet.length(),(double)_nInputs);

	// a likelihood is computed for every possible symbol vector
	tVector likelihoods(nSymbolVectors);

	tRange mPrecedentColumns(_startDetectionTime-_channelOrder+1,_startDetectionTime);
	tRange mMinus1PrecedentColumns(_startDetectionTime-_channelOrder+1,_startDetectionTime-1);

	// for each time instant
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
	{
		for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)
		{
			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			// the m-1 already detected symbol vectors are copied into the matrix:
			involvedSymbolVectors(rAllSymbolRows,mMinus1FirstColumns).inject(processedParticle->getSymbolVectors(mMinus1PrecedentColumns));

			for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.int2symbolsArray(iTestedVector,testedVector);

				// current tested vector is copied in the m-th position
				for(k=0;k<_nInputs;k++)
					involvedSymbolVectors(k,_channelOrder-1) = testedVector[k];

				likelihoods(iTestedVector) = processedParticle->getChannelMatrixEstimator(_estimatorIndex)->likelihood(observations.col(iObservationToBeProcessed),involvedSymbolVectors,noiseVariances[iObservationToBeProcessed]);
			} // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

			tVector probabilities(nSymbolVectors);
			try {
				// probabilities are computed by normalizing the likelihoods
				probabilities = Util::normalize(likelihoods);
			} catch (AllElementsNullException) {
				// if all the likelihoods are null
				probabilities = 1.0/(double)nSymbolVectors;
			}

			// one sample from the discrete distribution is taken
			iSampledVector = StatUtil::discrete_rnd(probabilities);

			// the above index is turned into a vector
			_alphabet.int2symbolsArray(iSampledVector,sampledVector);

			// sampled symbols are copied into the corresponding particle
			processedParticle->setSymbolVector(iObservationToBeProcessed,sampledVector);

			// channel matrix is estimated by means of the particle channel estimator
			processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),processedParticle->getSymbolVectors(mPrecedentColumns),noiseVariances[iObservationToBeProcessed]));

			processedParticle->setWeight(processedParticle->getWeight()* Util::sum(likelihoods));
		} // for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)

		mPrecedentColumns = mPrecedentColumns + 1;
		mMinus1PrecedentColumns = mMinus1PrecedentColumns + 1;

		_particleFilter->NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
// 			Resampling();
            _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);

	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

