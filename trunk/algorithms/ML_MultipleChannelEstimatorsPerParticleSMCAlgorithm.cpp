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

ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm::ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm)
{
// 	cout << "Probando..." << endl;
// 	for(int i=0;i<10;i++)
// 		cout << channelOrderPdf(1,0.5) << endl;
// 	exit(0);
}


ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm::~ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm()
{
}


void ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	int k,m,d,iSmoothingVector,nSmoothingVectors,Nm;
	int iSmoothingLag,iParticle,iSampledVector;
	int iSymbolVectorToBeProcessed;
// 	int m_previous;
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

// 		vector<vector<int> > indexes = _particleFilter.GetIndexesOfChannelOrders();
// 		for(int indexChannelOrder=0;indexChannelOrder<_candidateOrders.size();indexChannelOrder++)
// 		{
// 			cout << "De " << indexChannelOrder << "  hay " << indexes[indexChannelOrder].size() << endl;
// 		}

// 		iSymbolVectorToBeProcessed = _startDetectionSymbolVector+iObservationToBeProcessed-_startDetectionObservation;
//
// 			cout << "iObservationTo.. es " << iObservationToBeProcessed << endl;
		for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
		{
// 			cout << "iParticle es " << iParticle << endl;

			ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrder *> ( _particleFilter.GetParticle(iParticle));

// 			m_previous = processedParticle->GetChannelOrder();

			m = channelOrderPdf(processedParticle->GetChannelOrder(),P_CHANNEL_ORDER_PDF);

// 			cout << "Particula " << iParticle << " " << processedParticle->GetChannelOrder() << " -> " << m << endl;

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

// 			// channel matrix is estimated by means of the particle channel estimator
// 			processedParticle->SetChannelMatrix(iSymbolVectorToBeProcessed,(processedParticle->GetChannelMatrixEstimator())->NextMatrix(observations.col(iObservationToBeProcessed),processedParticle->GetSymbolVectors(iSymbolVectorToBeProcessed-m+1,iSymbolVectorToBeProcessed),noiseVariances[iObservationToBeProcessed]));

			// channel matrix is estimated with each of the channel estimators within the particle
			for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed,(processedParticle->GetChannelMatrixEstimator(iChannelOrder))->NextMatrix(observations.col(iObservationToBeProcessed),processedParticle->GetSymbolVectors(iObservationToBeProcessed-_candidateOrders[iChannelOrder]+1,iObservationToBeProcessed),noiseVariances[iObservationToBeProcessed]));
			}


			processedParticle->SetWeight(processedParticle->GetWeight()*Util::Sum(likelihoods));

// // 			cout << "Particula de orden " << processedParticle->GetChannelOrder() << "actualizada por " << Util::Sum(likelihoods) << endl;
//
// //             cout << "Terminada de procesar particula " << iParticle << endl;
//
		} // for(iParticle=0;iParticle<_nParticles;iParticle++)

// //         cout << "Instante " << iObservationToBeProcessed << endl;
//
// // 		NormalizeParticleGroups();
		_particleFilter.NormalizeWeights();

		cout << "Antes de resampling" << endl;

		vector<vector<int> > indexesPrueba = _particleFilter.GetIndexesOfChannelOrders();
		for(int indexChannelOrder=0;indexChannelOrder<_candidateOrders.size();indexChannelOrder++)
		{
			cout << "De " << indexChannelOrder << "  hay " << indexesPrueba[indexChannelOrder].size() << endl;
		}

		cout << "Los pesos" << endl << _particleFilter.GetWeightsVector() << endl;

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
		{
        	_resamplingAlgorithm->Resample(&_particleFilter);
		}

		cout << "Despues de resampling" << endl;

		indexesPrueba = _particleFilter.GetIndexesOfChannelOrders();
		for(int indexChannelOrder=0;indexChannelOrder<_candidateOrders.size();indexChannelOrder++)
		{
			cout << "De " << indexChannelOrder << "  hay " << indexesPrueba[indexChannelOrder].size() << endl;
		}

// 		for(int indexChannelOrder=0;indexChannelOrder<_candidateOrders.size();indexChannelOrder++)
// 		{
// 			cout << "De " << indexChannelOrder << "  hay " << _particleFilter.NparticlesOfChannelOrderIndex(indexChannelOrder) << endl;
// 		}

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
