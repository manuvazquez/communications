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
#include "ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm.h"

ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm::ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_particleFilter(nParticles,_candidateOrders)
{
}

void ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm::InitializeParticles()
{
    int iParticlePresentOrder,nParticlesPresentOrder,iParticle=0;
    int iChannelMatrixEstimator;

    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
    {
    // If the number of particles is not divisible by the number of channel orders, the remaining particles are assigned to the first channel orders (the + (......) term)
        nParticlesPresentOrder = _particleFilter.Nparticles()/_candidateOrders.size() + (_particleFilter.Nparticles() % _candidateOrders.size() > iChannelOrder);

        for(iParticlePresentOrder=0;iParticlePresentOrder<nParticlesPresentOrder;iParticlePresentOrder++,iParticle++)
        {
            // a clone of each of the channel matrix estimators is constructed...
            vector< ChannelMatrixEstimator * > thisParticleChannelMatrixEstimators(_candidateOrders.size());
            for(iChannelMatrixEstimator=0;iChannelMatrixEstimator<_candidateOrders.size();iChannelMatrixEstimator++)
                thisParticleChannelMatrixEstimators[iChannelMatrixEstimator] = _channelEstimators[iChannelMatrixEstimator]->Clone();

            // ... and passed within a vector to each particle
            _particleFilter.SetParticle(new ParticleWithChannelEstimationAndChannelOrder(1.0/(double)_particleFilter.Nparticles(),_N,_K,thisParticleChannelMatrixEstimators,_candidateOrders[iChannelOrder]),iParticle);
        }
    }
}

ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm::~ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm()
{
}


void ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
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

	// for each time instant
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{

		for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
		{

			ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrder *> ( _particleFilter.GetParticle(iParticle));

			m = channelOrderPdf(processedParticle->GetChannelOrder(),P_CHANNEL_ORDER_PDF);

			processedParticle->SetChannelOrder(m);

			// channel order dependent variables
			Nm = _N*m;
			d = m-1;
			nSmoothingVectors = (int) pow((double)_alphabet.Length(),(double)(_N*d));
			vector<tSymbol> testedSmoothingVector(_N*d);
			// it includes all symbol vectors involved in the smoothing
			tMatrix smoothingSymbolVectors(_N,m+d);

			// the m-1 already detected symbol vectors are copied into the matrix (just needed if m>1):
			if(m>1)
				smoothingSymbolVectors(allSymbolRows,tRange(0,m-2)).inject(processedParticle->GetSymbolVectors(iObservationToBeProcessed-m+1,iObservationToBeProcessed-1));

// 			cout << "Antes de iTestedVector" << endl;

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
					channelEstimatorClone = dynamic_cast < KalmanEstimator * > (processedParticle->GetChannelMatrixEstimator(_particleFilter.IndexFromChannelOrder(m))->Clone());

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
//
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
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledVector);

			// channel matrix is estimated with each of the channel estimators within the particle
			for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed,(processedParticle->GetChannelMatrixEstimator(iChannelOrder))->NextMatrix(observations.col(iObservationToBeProcessed),processedParticle->GetSymbolVectors(iObservationToBeProcessed-_candidateOrders[iChannelOrder]+1,iObservationToBeProcessed),noiseVariances[iObservationToBeProcessed]));
			}


			processedParticle->SetWeight(processedParticle->GetWeight()*Util::Sum(likelihoods));

		} // for(iParticle=0;iParticle<_nParticles;iParticle++)

		_particleFilter.NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
		{
        	_resamplingAlgorithm->Resample(&_particleFilter);
		}

	} // for each time instant
}

int ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm::channelOrderPdf(const int& m,const double& p)
{
	int alphabet[] = {-1,0,1};
	tVector probabilities(3);

	if(m==_minOrder)
	{
		probabilities(0) = 0;
		probabilities(1) = 1-p;
		probabilities(2) = p;
	}else if(m==_maxOrder) {
		probabilities(0) = p;
		probabilities(1) = 1-p;
		probabilities(2) = 0;		
	}else {
		probabilities(0) = p/2;
		probabilities(1) = 1-p;
		probabilities(2) = p/2;	
	}

	return m + alphabet[StatUtil::Discrete_rnd(1,probabilities)[0]];
}
