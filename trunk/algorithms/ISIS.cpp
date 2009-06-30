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
#include "ISIS.h"

// #define DEBUG3

ISIS::ISIS(string name, Alphabet alphabet, int L, int N, int frameLength, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, frameLength, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_particleFilter(nParticles)
{
}

void ISIS::InitializeParticles()
{
    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter.Capacity();iParticle++)
    {
		// a clone of each of the channel matrix estimators is constructed...
		vector< ChannelMatrixEstimator * > thisParticleChannelMatrixEstimators(_candidateOrders.size());
		for(uint iChannelMatrixEstimator=0;iChannelMatrixEstimator<_candidateOrders.size();iChannelMatrixEstimator++)
        {
			thisParticleChannelMatrixEstimators[iChannelMatrixEstimator] = _channelEstimators[iChannelMatrixEstimator]->Clone();
//             thisParticleChannelMatrixEstimators[iChannelMatrixEstimator]->setFirstEstimatedChannelMatrix(Util::ToMatrix(StatUtil::RandMatrix(_channelMeanVectors[iChannelMatrixEstimator],_channelCovariances[iChannelMatrixEstimator]),rowwise,_L));
        }

		// ... and passed within a vector to each particle
		_particleFilter.AddParticle(new ParticleWithChannelEstimationAndChannelOrderAPP(1.0/(double)_particleFilter.Capacity(),_N,_K,thisParticleChannelMatrixEstimators));
    }
}

void ISIS::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	int m,d,iSmoothingVector,nSmoothingVectors,Nm;
	int iSmoothingLag,iParticle,iSampledVector;
	uint iChannelOrder,k;
	vector<tSymbol> testedVector(_N),sampledVector(_N);
	double auxLikelihoodsProd,channelOrderAPPsNormConstant/*,newChannelOrderAPP*/;
	KalmanEstimator *auxChannelEstimator;
// 	double channelOrderAprioriProbability = 1.0/(double)_candidateOrders.size();

//     tVector newChannelOrderAPPs(_candidateOrders.size());
    double *newChannelOrderAPPs = new double[_candidateOrders.size()];

	// it selects all rows in the symbols Matrix
	tRange rAllSymbolRows(0,_N-1);

	int nSymbolVectors = (int) pow((double)_alphabet.length(),(double)_N);

	// a likelihood is computed for every possible symbol vector
	tVector likelihoods(nSymbolVectors);

	// for each time instant
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
#ifdef DEBUG2
        cout << "iObservationToBeProcessed = " << iObservationToBeProcessed << endl;
#endif

		for(iParticle=0;iParticle<_particleFilter.Capacity();iParticle++)
		{

            #ifdef DEBUG2
            cout << " ---------------- Particula " << iParticle << " ------------------------" << endl;
            #endif

			ParticleWithChannelEstimationAndChannelOrderAPP *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> ( _particleFilter.GetParticle(iParticle));

			for(int iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{

				likelihoods(iTestedVector) = 0.0;

				// the corresponding testing vector is generated from the index
				_alphabet.int2symbolsArray(iTestedVector,testedVector);

				for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
				{
					m = _candidateOrders[iChannelOrder];

					// channel order dependent variables
					Nm = _N*m;
// 					d = m-1;
					d = _maxOrder-1;
					nSmoothingVectors = (int) pow((double)_alphabet.length(),(double)(_N*d));
					vector<tSymbol> testedSmoothingVector(_N*d);
					// it includes all symbol vectors involved in the smoothing
					tMatrix smoothingSymbolVectors(_N,m+d);

					// the m-1 already detected symbol vectors are copied into the matrix (just needed if m>1):
					if(m>1)
						smoothingSymbolVectors(rAllSymbolRows,tRange(0,m-2)).inject(processedParticle->GetSymbolVectors(iObservationToBeProcessed-m+1,iObservationToBeProcessed-1));

					// current tested vector is copied in the m-th position
					for(k=0;k<_N;k++)
						smoothingSymbolVectors(k,m-1) = testedVector[k];

					// every possible smoothing sequence is tested
					for(iSmoothingVector=0;iSmoothingVector<nSmoothingVectors;iSmoothingVector++)
					{
						// a testing smoothing vector is generated for the index
						_alphabet.int2symbolsArray(iSmoothingVector,testedSmoothingVector);

						// symbols used for smoothing are copied into "smoothingSymbolVectors"
						for(k=0;k<testedSmoothingVector.size();k++)
							smoothingSymbolVectors((Nm+k)%_N,(Nm+k)/_N) = testedSmoothingVector[k];

						auxLikelihoodsProd = 1.0;

						// a clone of the channel estimator is generated because this must not be modified
						auxChannelEstimator = dynamic_cast < KalmanEstimator * > (processedParticle->GetChannelMatrixEstimator(iChannelOrder)->Clone());

						for(iSmoothingLag=0;iSmoothingLag<=d;iSmoothingLag++)
						{
							tRange mColumns(iSmoothingLag,iSmoothingLag+m-1);

							// the likelihood is computed and accumulated
							auxLikelihoodsProd *= auxChannelEstimator->likelihood(observations.col(iObservationToBeProcessed+iSmoothingLag),smoothingSymbolVectors(rAllSymbolRows,mColumns),noiseVariances[iObservationToBeProcessed+iSmoothingLag]);

							// a step in the Kalman Filter
							auxChannelEstimator->nextMatrix(observations.col(iObservationToBeProcessed+iSmoothingLag),smoothingSymbolVectors(rAllSymbolRows,mColumns),noiseVariances[iObservationToBeProcessed+iSmoothingLag]);
						} // for(iSmoothingLag=0;iSmoothingLag<=_d;iSmoothingLag++)

						// memory of the clone is freed
						delete auxChannelEstimator;

						// the likelihood is accumulated in the proper position of vector:
						likelihoods(iTestedVector) += auxLikelihoodsProd*processedParticle->GetChannelOrderAPP(iChannelOrder);

					} // for(iSmoothingVector=0;iSmoothingVector<nSmoothingVectors;iSmoothingVector++)

				} // for(iChannelOrder=0;iChannelOrder<_candidateOrders.size())

			} // for(int iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

			tVector probabilities(nSymbolVectors);
			try {
				// probabilities are computed by normalizing the likelihoods
				probabilities = Util::Normalize(likelihoods);
			} catch (AllElementsNullException) {
				// if all the likelihoods are null
				probabilities = 1.0/(double)nSymbolVectors;
			}

            #ifdef DEBUG
                cout << "Likelihoods" << endl << likelihoods << endl << "Probabilidades" << endl << probabilities << endl;
                cout << "Simbolos transmitidos" << endl << _simbolos.col(iObservationToBeProcessed) << endl;
            #endif

			// one sample from the discrete distribution is taken
			iSampledVector = StatUtil::discrete_rnd(probabilities);

			// the above index is turned into a vector
			_alphabet.int2symbolsArray(iSampledVector,sampledVector);

			// sampled symbols are copied into the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledVector);

            #ifdef DEBUG
                for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
                {
                    cout << "memoria " << _candidateOrders[iChannelOrder] << ": " << processedParticle->GetChannelOrderAPP(iChannelOrder) << " ";
                }
                cout << endl;
            #endif

            channelOrderAPPsNormConstant = 0.0;

			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				auxChannelEstimator = dynamic_cast < KalmanEstimator * > (processedParticle->GetChannelMatrixEstimator(iChannelOrder));

                // for efficiency's sake
                tMatrix involvedSymbolVectors = processedParticle->GetSymbolVectors (iObservationToBeProcessed-_candidateOrders[iChannelOrder]+1,iObservationToBeProcessed);

                // the a posteriori probability for each channel order must be updated with the previous app for this order and the likelihood at the present instant with the sampled symbol vector
                newChannelOrderAPPs[iChannelOrder] = processedParticle->GetChannelOrderAPP(iChannelOrder)*
                auxChannelEstimator->likelihood(observations.col(iObservationToBeProcessed),involvedSymbolVectors,noiseVariances[iObservationToBeProcessed]);

                channelOrderAPPsNormConstant += newChannelOrderAPPs[iChannelOrder];

// 				processedParticle->SetChannelOrderAPP(newChannelOrderAPP,iChannelOrder);

				// channel matrix is estimated with each of the channel estimators within the particle
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed,auxChannelEstimator->nextMatrix(observations.col(iObservationToBeProcessed),involvedSymbolVectors,noiseVariances[iObservationToBeProcessed]));
			}

            if(channelOrderAPPsNormConstant!=0)
                // all the channel order a posteriori probabilities are normalized by the previously computed constant
                for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
                    processedParticle->SetChannelOrderAPP(newChannelOrderAPPs[iChannelOrder]/channelOrderAPPsNormConstant,iChannelOrder);


            #ifdef DEBUG
                cout << "Despues" << endl;
                for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
                {
                    cout << "memoria " << _candidateOrders[iChannelOrder] << ": " << processedParticle->GetChannelOrderAPP(iChannelOrder) << " ";
                }
                cout << endl;
            #endif

			processedParticle->SetWeight(processedParticle->GetWeight()*Util::Sum(likelihoods));

		} // for(iParticle=0;iParticle<_nParticles;iParticle++)

		_particleFilter.NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
		{
        	_resamplingAlgorithm->ResampleWhenNecessary(&_particleFilter);
		}

	} // for each time instant

	delete[] newChannelOrderAPPs;
}

vector<tMatrix> ISIS::GetEstimatedChannelMatrices()
{
	return vector<tMatrix>(0);
}
