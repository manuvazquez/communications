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

ML_UnknownChannelOrderSMCAlgorithm::ML_UnknownChannelOrderSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles,ResamplingAlgorithm *resamplingAlgorithm,ResamplingAlgorithm *resamplingAlgorithm2,tMatrix simbolosVerdaderos): UnknownChannelOrderSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles,resamplingAlgorithm),_simbolosVerdaderos(simbolosVerdaderos),_resamplingAlgorithm2(resamplingAlgorithm2)
{
}

void ML_UnknownChannelOrderSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
//     cout << "El tamaño de las observaciones es " << observations.cols() << endl;
	int k,m,d,iSmoothingVector,nSmoothingVectors,Nm;
	int iSmoothingLag,iParticle,iSampledVector;
	int iSymbolVectorToBeProcessed;
	vector<tSymbol> testedVector(_N),sampledVector(_N);
	double auxLikelihoodsProd;
	KalmanEstimator *channelEstimatorClone;

	// it selects all rows in the symbols Matrix
	tRange allSymbolRows(0,_N-1);

	int nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);

	// a likelihood is computed for every possible symbol vector
	tVector likelihoods(nSymbolVectors);

    cout << "_K es " << _K << endl;

	// for each time instant
	for(int iObservationToBeProcessed=_startDetectionObservation;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{

		iSymbolVectorToBeProcessed = _startDetectionSymbolVector+iObservationToBeProcessed-_startDetectionObservation;

// 		cout << "Observacion procesada: " << iObservationToBeProcessed << endl;

		for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
		{
			ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrder *> ( _particleFilter.GetParticle(iParticle));

			m = processedParticle->GetChannelOrder();

//             cout << "El orden de la particula es " << m << endl;

			// channel order dependent variables
			Nm = _N*m;
			d = m-1;
			nSmoothingVectors = (int) pow((double)_alphabet.Length(),(double)(_N*d));
			vector<tSymbol> testedSmoothingVector(_N*d);
			// it includes all symbol vectors involved in the smoothing
			tMatrix smoothingSymbolVectors(_N,m+d);

			// the m-1 already detected symbol vectors are copied into the matrix:
			smoothingSymbolVectors(allSymbolRows,tRange(0,m-2)).inject(processedParticle->GetSymbolVectors(iSymbolVectorToBeProcessed-m+1,iSymbolVectorToBeProcessed-1));

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
			processedParticle->SetSymbolVector(iSymbolVectorToBeProcessed,sampledVector);

			// channel matrix is estimated by means of the particle channel estimator
			processedParticle->SetChannelMatrix(iSymbolVectorToBeProcessed,(processedParticle->GetChannelMatrixEstimator())->NextMatrix(observations.col(iObservationToBeProcessed),processedParticle->GetSymbolVectors(iSymbolVectorToBeProcessed-m+1,iSymbolVectorToBeProcessed),noiseVariances[iObservationToBeProcessed]));

			processedParticle->SetWeight(processedParticle->GetWeight()* Util::Sum(likelihoods));

// 			cout << "Particula de orden " << processedParticle->GetChannelOrder() << "actualizada por " << Util::Sum(likelihoods) << endl;

//             cout << "Terminada de procesar particula " << iParticle << endl;

		} // for(iParticle=0;iParticle<_nParticles;iParticle++)

//         cout << "Instante " << iObservationToBeProcessed << endl;

// 		NormalizeParticleGroups();
		_particleFilter.NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
		{
			if(iObservationToBeProcessed<15)
            {
// 				ResamplingByParticleGroups();
                _resamplingAlgorithm->Resample(&_particleFilter);
            }
			else
                _resamplingAlgorithm2->Resample(&_particleFilter);
// 				Resampling();
		}

	} // for each time instant
}

