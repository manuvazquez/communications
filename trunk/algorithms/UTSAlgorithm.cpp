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
#include "UTSAlgorithm.h"

// #define DEBUG

UTSAlgorithm::UTSAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, double ARcoefficient, double samplingVariance, double ARprocessVariance): UPSPBasedSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm, ARcoefficient, samplingVariance, ARprocessVariance)
{
}

void UTSAlgorithm::InitializeParticles()
{
	vector<ChannelMatrixEstimator *> channelEstimatorsClone(_channelEstimators.size());
	for(uint i=0;i<_candidateOrders.size();i++)
		channelEstimatorsClone[i] = _channelEstimators[i]->Clone();

	// we begin with only one particle
	ParticleWithChannelEstimationAndChannelOrderAPP *particle = new ParticleWithChannelEstimationAndChannelOrderAPP(1.0,_N,_K+_d,channelEstimatorsClone);

	particle->SetSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);

	// the available APP's just before the _startDetectionTime instant are copied into the particle
	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
		particle->SetChannelOrderAPP(_channelOrderAPPs(iChannelOrder,_startDetectionTime-1),iChannelOrder);

	_particleFilter->AddParticle(particle);
}

void UTSAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	uint nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);
	tRange rMaxChannelOrderMinus1FirstColumns(0,_maxOrder-2),rAll;
	vector<tSymbol> testedVector(_N);
	tVector computedObservations(_L);
	vector<tVector> computedObservationsVector(_candidateOrders.size());
	int iCandidate,m,iBestUnnormalizedChannelOrderAPP,k,iParticle,iTestedVector;
	uint iChannelOrder;
	ParticleWithChannelEstimationAndChannelOrderAPP *processedParticle;
	double likelihood;

	typedef struct{
		int fromParticle;
		tMatrix symbolVectorsMatrix;
        int iBestChannelOrder;
        vector<tVector> computedObservationsVector;
		tVector unnormalizedChannelOrderAPPs;
		double likelihood;
		double weight;
	}tParticleCandidate;

	tParticleCandidate *particleCandidates = new tParticleCandidate[_particleFilter->Capacity()*nSymbolVectors];

    // "symbolVectorsMatrix" will contain all the symbols involved in the current observation
    tMatrix symbolVectorsMatrix(_N,_maxOrder);
    tVector symbolsVector(_NmaxOrder);

    int lastSymbolVectorStart = _NmaxOrder - _N;

	tRange rMaxChannelOrderMinus1PrecedentColumns(_startDetectionTime-_maxOrder+1,_startDetectionTime-1);

	vector<bool> activeCandidateOrders(_candidateOrders.size(),true);
	int iBestChannelOrder = 0,timesBestChannelOrder = 0;

	// at first, there is only one particle
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K+_d;iObservationToBeProcessed++)
	{
#ifdef DEBUG
// 		cout << iObservationToBeProcessed << endl;
#endif

		// it keeps track of the place where a new tParticleCandidate will be stored within the array
		iCandidate = 0;

		// the candidates from all the particles are generated
		for(iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
		{
			processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->GetParticle(iParticle));

			symbolVectorsMatrix(rAll,rMaxChannelOrderMinus1FirstColumns).inject(processedParticle->GetSymbolVectors(rMaxChannelOrderMinus1PrecedentColumns));
			symbolsVector = Util::ToVector(symbolVectorsMatrix,columnwise);

			for(iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.IntToSymbolsArray(iTestedVector,testedVector);

				// current tested vector is copied in the m-th position
				for(k=0;k<_N;k++)
					symbolVectorsMatrix(k,_maxOrder-1) = symbolsVector(lastSymbolVectorStart+k) = testedVector[k];

				likelihood = 0.0;

				particleCandidates[iCandidate].unnormalizedChannelOrderAPPs = tVector(_candidateOrders.size());

				for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
				{
					m = _candidateOrders[iChannelOrder];

					tMatrix estimatedChannelMatrix = processedParticle->GetChannelMatrixEstimator(iChannelOrder)->LastEstimatedChannelMatrix();
					estimatedChannelMatrix *= _ARcoefficient;

					// computedObservations = estimatedChannelMatrix * symbolVectorsMatrix(:)
					Blas_Mat_Vec_Mult(estimatedChannelMatrix,symbolsVector(tRange(_NmaxOrder-m*_N,_NmaxOrder-1)),computedObservations);

					particleCandidates[iCandidate].unnormalizedChannelOrderAPPs(iChannelOrder) = processedParticle->GetChannelOrderAPP(iChannelOrder)*StatUtil::NormalPdf(observations.col(iObservationToBeProcessed),computedObservations,noiseVariances[iObservationToBeProcessed]);

					likelihood += particleCandidates[iCandidate].unnormalizedChannelOrderAPPs(iChannelOrder);

					computedObservationsVector[iChannelOrder] = computedObservations;
				}

				// if the likelihood is zero, we don't generate a candidate for this particle and this symbol vector
				if(likelihood==0.0)
					continue;
// 					throw RuntimeException("UTSAlgorithm::Process: likelihood is zero.");

				Util::Max(particleCandidates[iCandidate].unnormalizedChannelOrderAPPs,iBestUnnormalizedChannelOrderAPP);

				particleCandidates[iCandidate].fromParticle = iParticle;
				particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;
                particleCandidates[iCandidate].iBestChannelOrder = iBestUnnormalizedChannelOrderAPP;
                particleCandidates[iCandidate].computedObservationsVector = computedObservationsVector;
                particleCandidates[iCandidate].likelihood = likelihood;
				particleCandidates[iCandidate].weight = processedParticle->GetWeight()*likelihood;

				iCandidate++;
			} // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

		} // for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)

		// a vector of size the number of generated candidates is declared...
		tVector weights(iCandidate);

		// ...to store their weights
		for(int i=0;i<iCandidate;i++)
			weights(i) = particleCandidates[i].weight;

		// the candidates that are going to give particles are selected
		vector<int> indexesSelectedCandidates = _resamplingAlgorithm->ObtainIndexes(_particleFilter->Capacity(),weights);

		// every survivor candidate is associated with an old particle
		vector<int> indexesParticles(indexesSelectedCandidates.size());
		for(uint i=0;i<indexesSelectedCandidates.size();i++)
        {
			indexesParticles[i] = particleCandidates[indexesSelectedCandidates[i]].fromParticle;
            _particlesBestChannelOrders[i] = particleCandidates[indexesSelectedCandidates[i]].iBestChannelOrder;
        }

		// the chosen particles are kept without modification (yet)
		_particleFilter->KeepParticles(indexesParticles);

		// every surviving particle is modified according to what it says its corresponding candidate
		for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
		{
			processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->GetParticle(iParticle));

			// sampled symbols are copied into the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.col(_maxOrder-1));

			for(uint iChannelOrder=0;iChannelOrder<processedParticle->NchannelMatrixEstimators();iChannelOrder++)
			{
				// channel matrix is estimated by means of the particle channel estimator
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed,processedParticle->GetChannelMatrixEstimator(iChannelOrder)->NextMatrix(observations.col(iObservationToBeProcessed),particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix(rAll,tRange(_maxOrder-_candidateOrders[iChannelOrder],_maxOrder-1)),noiseVariances[iObservationToBeProcessed]));

				processedParticle->SetChannelOrderAPP(particleCandidates[indexesSelectedCandidates[iParticle]].unnormalizedChannelOrderAPPs(iChannelOrder)/particleCandidates[indexesSelectedCandidates[iParticle]].likelihood,iChannelOrder);
			}

			processedParticle->SetWeight(particleCandidates[indexesSelectedCandidates[iParticle]].weight);

		} // for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)

		_particleFilter->NormalizeWeights();

		processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->GetBestParticle());

        for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
	        _channelOrderAPPs(iChannelOrder,iObservationToBeProcessed) = processedParticle->GetChannelOrderAPP(iChannelOrder);

		if(_particlesBestChannelOrders[_particleFilter->iBestParticle()]==iBestChannelOrder)
			timesBestChannelOrder++;
		else
		{
			iBestChannelOrder = _particlesBestChannelOrders[_particleFilter->iBestParticle()];
			timesBestChannelOrder = 0;
		}

#ifdef DEBUG
		if(timesBestChannelOrder > 10)
			cout << "iBestChannelOrder = " << iBestChannelOrder << " ha sido el mejor durante " << timesBestChannelOrder << endl;
#endif

#ifdef DEBUG
		cout << _channelOrderAPPs.col(iObservationToBeProcessed);
		cout << "---------------------" << endl;
#endif

		rMaxChannelOrderMinus1PrecedentColumns = rMaxChannelOrderMinus1PrecedentColumns + 1;

	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K+_d;iObservationToBeProcessed++)

#ifdef DEBUG
// 	cout << _channelOrderAPPs(rAll,tRange(250,260)) << endl;
#endif

	delete[] particleCandidates;
}

