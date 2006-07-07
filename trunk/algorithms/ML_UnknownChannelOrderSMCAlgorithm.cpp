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
#include "ML_UnknownChannelOrderSMCAlgorithm.h"

ML_UnknownChannelOrderSMCAlgorithm::ML_UnknownChannelOrderSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int firstObservationIndex, int smoothingLag, int nParticles, ResamplingCriterion resamplingCriterion, StdResamplingAlgorithm resamplingAlgorithm): UnknownChannelOrderSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, firstObservationIndex, smoothingLag, nParticles, resamplingCriterion, resamplingAlgorithm)
{
}


ML_UnknownChannelOrderSMCAlgorithm::~ML_UnknownChannelOrderSMCAlgorithm()
{
}


void ML_UnknownChannelOrderSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	cout << "En process..." << endl;
	int k,m,d,iSmoothingVector,nSmoothingVectors,Nm;
	int iSmoothingLag,iParticle,iSampledVector;
	vector<tSymbol> testedVector(_N),sampledVector(_N);
	double auxLikelihoodsProd;
	KalmanEstimator *channelEstimatorClone;
// 	tRange mFirstColumns(0,_m-1);

	// it selects all rows in the symbols Matrix
	tRange allSymbolRows(0,_N-1);

	// it includes all symbol vectors involved in the smoothing
// 	tMatrix smoothingSymbolVectors(_N,_m+_d);

	int nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);
// 	int nSmoothingVectors = (int) pow((double)_alphabet.Length(),(double)(_N*_d));
	
	// a likelihood is computed for every possible symbol vector
	tVector likelihoods(nSymbolVectors);

	// for each time instant
	for(int iObservationToBeProcessed=_startDetectionObservation;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
		cout << "Observacion procesada: " << iObservationToBeProcessed << endl;

// 		tRange mPrecedentColumns(iObservationToBeProcessed-_m+1,iObservationToBeProcessed);
		for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
		{
			ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrder *> ( _particleFilter.GetParticle(iParticle));

			m = processedParticle->GetChannelOrder();

			// channel order dependent variables
			Nm = _N*m;
			d = m-1;
			nSmoothingVectors = (int) pow((double)_alphabet.Length(),(double)(_N*d));
			vector<tSymbol> testedSmoothingVector(_N*d);
			// it includes all symbol vectors involved in the smoothing	
			tMatrix smoothingSymbolVectors(_N,m+d);

			// the m-1 already detected symbol vectors are copied into the matrix:
			smoothingSymbolVectors(allSymbolRows,tRange(0,m-2)).inject(processedParticle->GetSymbolVectors(_startDetectionSymbolVector-_startDetectionObservation+iObservationToBeProcessed-m+1,_startDetectionSymbolVector-_startDetectionObservation+iObservationToBeProcessed-1));

			for(int iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.IntToSymbolsArray(iTestedVector,testedVector);
				
				// current tested vector is copied in the m-th position
				for(k=0;k<_N;k++)
					smoothingSymbolVectors(k,m-1) = testedVector[k];

				likelihoods(iTestedVector) = 0.0;

				// every possible smoothing sequence is tested
				for(iSmoothingVector=0;iSmoothingVector<nSmoothingVectors;iSmoothingVector++)
				{
					// a testing smoothing vector is generated for the index
					_alphabet.IntToSymbolsArray(iSmoothingVector,testedSmoothingVector);

					// symbols used for smoothing are copied into "smoothingSymbolVectors"
					for(k=0;k<testedSmoothingVector.size();k++)
						smoothingSymbolVectors((Nm+k)%_N,(Nm+k)/_N) = testedSmoothingVector[k];

					auxLikelihoodsProd = 1.0;

					// a clone of the channel estimator is generated because this must not be modified
					channelEstimatorClone = dynamic_cast < KalmanEstimator * > (processedParticle->GetChannelMatrixEstimator()->Clone());

					for(iSmoothingLag=0;iSmoothingLag<=d;iSmoothingLag++)
					{
						tRange mColumns(iSmoothingLag,iSmoothingLag+m-1);

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

			tVector probabilities(nSymbolVectors);
			try {
				// probabilities are computed by normalizing the likelihoods
				probabilities = Util::Normalize(likelihoods);
			} catch (AllElementsNullException) {
				// if all the likelihoods are null
				probabilities = 1.0/(double)nSymbolVectors;
			}

			// one sample from the discrete distribution is taken
			int iSampledVector = (StatUtil::Discrete_rnd(1,probabilities))[0];

			// the above index is turned into a vector
			_alphabet.IntToSymbolsArray(iSampledVector,sampledVector);

			// sampled symbols are copied into the corresponding particle
			processedParticle->SetSymbolVector(_startDetectionSymbolVector-_startDetectionObservation+iObservationToBeProcessed,sampledVector);
						
			// channel matrix is estimated by means of the particle channel estimator
			processedParticle->SetChannelMatrix(_startDetectionSymbolVector-_startDetectionObservation+iObservationToBeProcessed,(processedParticle->GetChannelMatrixEstimator())->NextMatrix(observations.col(iObservationToBeProcessed),processedParticle->GetSymbolVectors(_startDetectionSymbolVector-_startDetectionObservation+iObservationToBeProcessed-m+1,_startDetectionSymbolVector-_startDetectionObservation+iObservationToBeProcessed),noiseVariances[iObservationToBeProcessed]));

			processedParticle->SetWeight(processedParticle->GetWeight()* Util::Sum(likelihoods));

		} // for(iParticle=0;iParticle<_nParticles;iParticle++)

// 		cout << "Terminé" << endl;
// 		for(int holita=0;holita<_nCandidateOrders;holita++)
// 			cout << "orden " << _candidateOrders[holita] << " -> "  << "indice " << holita << endl;
// 		exit(0);

		NormalizeParticleGroups();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
			Resampling();
	} // for each time instant
	cout << "Terminé" << endl;
	exit(0);
}

