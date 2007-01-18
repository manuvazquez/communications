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
#include "ML_SMCAlgorithm.h"

// #define DEBUG

ML_SMCAlgorithm::ML_SMCAlgorithm(string name, Alphabet alphabet,int L,int N, int K,int m, ChannelMatrixEstimator *channelEstimator, tMatrix preamble, int smoothingLag, int nParticles,StdResamplingAlgorithm resamplingAlgorithm): SMCAlgorithm(name, alphabet, L, N, K,m,  channelEstimator, preamble, smoothingLag, nParticles,resamplingAlgorithm)
{
}

void ML_SMCAlgorithm::Process(const tMatrix &observations, vector< double > noiseVariances)
{
	int k,iSmoothingVector;
	int iSmoothingLag,iParticle,iSampledVector;
	vector<tSymbol> testedVector(_N),testedSmoothingVector(_N*_d),sampledVector(_N);
	double auxLikelihoodsProd;
	ChannelMatrixEstimator *channelEstimatorClone;
	tRange mFirstColumns(0,_m-1);

	// it selects all rows in the symbols Matrix
	tRange allSymbolRows(0,_N-1);

	// it includes all symbol vectors involved in the smoothing
	tMatrix smoothingSymbolVectors(_N,_m+_d);

	int nObservations = observations.cols();
	int nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);
	int nSmoothingVectors = (int) pow((double)_alphabet.Length(),(double)(_N*_d));

	// a likelihood is computed for every possible symbol vector
	tVector likelihoods(nSymbolVectors);

	// for each time instant
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
		#ifdef DEBUG
			cout << "Observacion procesada: " << iObservationToBeProcessed << endl;
		#endif

		tRange mPrecedentColumns(iObservationToBeProcessed-_m+1,iObservationToBeProcessed);
		for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
		{
			#ifdef DEBUG
				cout << "Particula: " << iParticle << endl;
			#endif

			ParticleWithChannelEstimation *processedParticle = _particleFilter.GetParticle(iParticle);

			// the m-1 already detected symbol vectors are copied into the matrix:
			smoothingSymbolVectors(allSymbolRows,mFirstColumns).inject(processedParticle->GetSymbolVectors(mPrecedentColumns));

			for(int iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.IntToSymbolsArray(iTestedVector,testedVector);

				// current tested vector is copied in the m-th position
				for(k=0;k<_N;k++)
					smoothingSymbolVectors(k,_m-1) = testedVector[k];

				likelihoods(iTestedVector) = 0.0;

				// every possible smoothing sequence is tested
				for(iSmoothingVector=0;iSmoothingVector<nSmoothingVectors;iSmoothingVector++)
				{
					// a testing smoothing vector is generated for the index
					_alphabet.IntToSymbolsArray(iSmoothingVector,testedSmoothingVector);

					// symbols used for smoothing are copied into "smoothingSymbolVectors"
					for(k=0;k<testedSmoothingVector.size();k++)
						smoothingSymbolVectors((_Nm+k)%_N,(_Nm+k)/_N) = testedSmoothingVector[k];

					auxLikelihoodsProd = 1.0;

					#ifdef DEBUG
						cout << "Before clonig the channel estimator." << endl;
					#endif

					// a clone of the channel estimator is generated because this must not be modified
					channelEstimatorClone = processedParticle->GetChannelMatrixEstimator()->Clone();

					#ifdef DEBUG
						cout << "After clonig the channel estimator." << endl;
					#endif

					for(iSmoothingLag=0;iSmoothingLag<=_d;iSmoothingLag++)
					{
						tRange mColumns(iSmoothingLag,iSmoothingLag+_m-1);

						// the likelihood is computed and accumulated
						auxLikelihoodsProd *= channelEstimatorClone->Likelihood(observations.col(iObservationToBeProcessed+iSmoothingLag),smoothingSymbolVectors(allSymbolRows,mColumns),noiseVariances[iObservationToBeProcessed+iSmoothingLag]);

						#ifdef DEBUG
							cout << "Despues de llamar a likelihood." << endl;
						#endif

						// a step in the Kalman Filter
						channelEstimatorClone->NextMatrix(observations.col(iObservationToBeProcessed+iSmoothingLag),smoothingSymbolVectors(allSymbolRows,mColumns),noiseVariances[iObservationToBeProcessed+iSmoothingLag]);
					} // for(iSmoothingLag=0;iSmoothingLag<=_d;iSmoothingLag++)

					// memory of the clone is freed
					delete channelEstimatorClone;

					// the likelihood is accumulated in the proper position of vector:
					likelihoods(iTestedVector) += auxLikelihoodsProd;
				} // for(iSmoothingVector=0;iSmoothingVector<nSmoothingVectors;iSmoothingVector++)

			} // for(int iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

			tVector probabilities(nSymbolVectors);
			try {
				// probabilities are computed by normalizing the likelihoods
				probabilities = Util::Normalize(likelihoods);
			} catch (AllElementsNullException) {
				// if all the likelihoods are null
				probabilities = 1.0/(double)nSymbolVectors;
			}

			// one sample from the discrete distribution is taken
			int iSampledVector = StatUtil::Discrete_rnd(probabilities);

			// the above index is turned into a vector
			_alphabet.IntToSymbolsArray(iSampledVector,sampledVector);

			// sampled symbols are copied into the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledVector);

			// channel matrix is estimated by means of the particle channel estimator
			processedParticle->SetChannelMatrix(iObservationToBeProcessed,(processedParticle->GetChannelMatrixEstimator())->NextMatrix(observations.col(iObservationToBeProcessed),processedParticle->GetSymbolVectors(mPrecedentColumns),noiseVariances[iObservationToBeProcessed]));

			processedParticle->SetWeight(processedParticle->GetWeight()* Util::Sum(likelihoods));

		} // for(iParticle=0;iParticle<_nParticles;iParticle++)

		_particleFilter.NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
// 			Resampling();
            _resamplingAlgorithm.Resample(&_particleFilter);

	} // for each time instant
}

