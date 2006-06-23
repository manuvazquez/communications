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

ML_SMCAlgorithm::ML_SMCAlgorithm(string name, Alphabet alphabet, int K, ChannelMatrixEstimator *channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingCriterion resamplingCriterion,StdResamplingAlgorithm resamplingAlgorithm): SMCAlgorithm(name, alphabet, K,  channelEstimator, preamble, smoothingLag, nParticles, resamplingCriterion,resamplingAlgorithm)
{
// 	const type_info &info = typeid(channelEstimator);
// 	cout << "El nombre del tipo es:" << typeid(channelEstimator).name() << endl;
// 	char c;
// 	cin >> c;
}


void ML_SMCAlgorithm::Process(const tMatrix &observations, vector< double > noiseVariances)
{
// 	cout << "En Process" << endl;

	int k,iSmoothingVector;
	int iSmoothingLag,iParticle,iSampledVector;
	vector<tSymbol> testedVector(_N),testedSmoothingVector(_N*_d),sampledVector(_N);
	double auxLikelihoodsProd;
	KalmanEstimator *channelEstimatorClone;
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
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_endDetectionTime;iObservationToBeProcessed++)
	{
// 		cout << "Observacion procesada: " << iObservationToBeProcessed << endl;

		tRange mPrecedentColumns(iObservationToBeProcessed-_m+1,iObservationToBeProcessed);
		for(iParticle=0;iParticle<_nParticles;iParticle++)
		{
			// the m-1 already detected symbol vectors are copied into the matrix:
			smoothingSymbolVectors(allSymbolRows,mFirstColumns).inject(_particles[iParticle]->GetSymbolVectors(mPrecedentColumns));

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

					// a clone of the channel estimator is generated because this must not be modified
					channelEstimatorClone = dynamic_cast < KalmanEstimator * > (_particles[iParticle]->GetChannelMatrixEstimator()->Clone());

					for(iSmoothingLag=0;iSmoothingLag<=_d;iSmoothingLag++)
					{
						tRange mColumns(iSmoothingLag,iSmoothingLag+_m-1);

						// the likelihood is computed and accumulated
						auxLikelihoodsProd *= channelEstimatorClone->Likelihood(observations.col(iObservationToBeProcessed+iSmoothingLag),smoothingSymbolVectors(allSymbolRows,mColumns),noiseVariances[iObservationToBeProcessed+iSmoothingLag]);

						// a step in the Kalman Filter
						channelEstimatorClone->NextMatrix(observations.col(iObservationToBeProcessed+iSmoothingLag),smoothingSymbolVectors(allSymbolRows,mColumns),noiseVariances[iObservationToBeProcessed+iSmoothingLag]);
					} // for(iSmoothingLag=0;iSmoothingLag<=_d;iSmoothingLag++)

					// memory of the clone is freed
					delete channelEstimatorClone;

					// the likelihood is accumulated in the proper position of vector:
					likelihoods(iTestedVector) += auxLikelihoodsProd;
				} // for(iSmoothingVector=0;iSmoothingVector<nSmoothingVectors;iSmoothingVector++)

			} // for(int iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

			// probabilities are computed by normalizing the likelihoods
			tVector probabilities = Util::Normalize(likelihoods);

			// one sample from the discrete distribution is taken
			int iSampledVector = (StatUtil::Discrete_rnd(1,probabilities))[0];

			// the above index is turned into a vector
			_alphabet.IntToSymbolsArray(iSampledVector,sampledVector);

			// sampled symbols are copied into the corresponding particle
			_particles[iParticle]->SetSymbolVector(iObservationToBeProcessed,sampledVector);
						
			// channel matrix is estimated by means of the particle channel estimator
			_particles[iParticle]->SetChannelMatrix(iObservationToBeProcessed,(_particles[iParticle]->GetChannelMatrixEstimator())->NextMatrix(observations.col(iObservationToBeProcessed),_particles[iParticle]->GetSymbolVectors(mPrecedentColumns),noiseVariances[iObservationToBeProcessed]));

			_particles[iParticle]->SetWeight(_particles[iParticle]->GetWeight()* Util::Sum(likelihoods));

		} // for(iParticle=0;iParticle<_nParticles;iParticle++)

		NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_endDetectionTime-1))
			this->Resampling(iObservationToBeProcessed);

	} // for each time instant
}

