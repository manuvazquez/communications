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
#include "UPSPBasedSMCAlgorithm.h"

// #define DEBUG

UPSPBasedSMCAlgorithm::UPSPBasedSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_particleFilter(new ParticleFilter(nParticles)),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_particlesBestChannelOrders(nParticles)
{
}


UPSPBasedSMCAlgorithm::~UPSPBasedSMCAlgorithm()
{
	#ifdef DEBUG2
		cout << "holita" << endl;
	#endif
	delete _particleFilter;

	#ifdef DEBUG2
		cout << "adios" << endl;
	#endif
}

void UPSPBasedSMCAlgorithm::InitializeParticles()
{
	#ifdef DEBUG
		cout <<  __FILE__  << "(line " << __LINE__ << ")" << endl;
		cout << "Al principio de InitializeParticles" << endl;
	#endif

	vector<ChannelMatrixEstimator *> channelEstimatorsClone(_channelEstimators.size());
	for(uint i=0;i<_candidateOrders.size();i++)
		channelEstimatorsClone[i] = _channelEstimators[i]->Clone();

	// we begin with only one particle
	_particleFilter->AddParticle(new ParticleWithChannelEstimation(1.0,_N,_K+_d,channelEstimatorsClone));
	_particleFilter->GetParticle(0)->SetSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);

	#ifdef DEBUG
		cout <<  __FILE__  << "(line " << __LINE__ << ")" << endl;
		cout << "Al final de InitializeParticles" << endl;
	#endif
}

void UPSPBasedSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	#ifdef DEBUG2
		cout << "Al principio de Process" << endl;
	#endif
	uint nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);
	tRange rMaxChannelOrderMinus1FirstColumns(0,_maxOrder-2),rAllSymbolRows(0,_N-1);
	vector<tSymbol> testedVector(_N);
	tVector computedObservations(_L),error(_L);
	vector<double> costs(_candidateOrders.size());
	vector<tVector> computedObservationsVector(_candidateOrders.size());

	#ifdef DEBUG
		cout << "Despues de inicializar las variables de process." << endl;
	#endif


	typedef struct{
		int fromParticle;
		tMatrix symbolVectorsMatrix;
        int iBestChannelOrder;
		double weight;
	}tParticleCandidate;

	tParticleCandidate *particleCandidates = new tParticleCandidate[_particleFilter->Capacity()*nSymbolVectors] ;

    // "symbolVectorsMatrix" will contain all the symbols involved in the current observation
    tMatrix symbolVectorsMatrix(_N,_maxOrder);
    tVector symbolsVector(_NmaxOrder);

    int lastSymbolVectorStart = _NmaxOrder - _N;

	// at first, there is only one particle
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K+_d;iObservationToBeProcessed++)
	{
		tRange rMaxChannelOrderMinus1PrecedentColumns(iObservationToBeProcessed-_maxOrder+1,iObservationToBeProcessed-1);

		// it keeps track of the place where a new tParticleCandidate will be stored within the array
		int iCandidate = 0;

		// the candidates from all the particles are generated
		for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
		{
			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			symbolVectorsMatrix(rAllSymbolRows,rMaxChannelOrderMinus1FirstColumns).inject(processedParticle->GetSymbolVectors(rMaxChannelOrderMinus1PrecedentColumns));
			symbolsVector = Util::ToVector(symbolVectorsMatrix,columnwise);

			for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.IntToSymbolsArray(iTestedVector,testedVector);

				// current tested vector is copied in the m-th position
				for(int k=0;k<_N;k++)
					symbolVectorsMatrix(k,_maxOrder-1) = symbolsVector(lastSymbolVectorStart+k) = testedVector[k];


				for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
				{
					int m = _candidateOrders[iChannelOrder];

					tMatrix estimatedChannelMatrix = processedParticle->GetChannelMatrixEstimator(iChannelOrder)->LastEstimatedChannelMatrix();
					estimatedChannelMatrix *= _ARcoefficient;

					// computedObservations = estimatedChannelMatrix * symbolVectorsMatrix(:)
					Blas_Mat_Vec_Mult(estimatedChannelMatrix,symbolsVector(tRange(_NmaxOrder-m*_N,_NmaxOrder-1)),computedObservations);

					// error = observations.col(iObservationToBeProcessed) - computedObservations
					Util::Add(tVector(observations.col(iObservationToBeProcessed)),tVector(computedObservations),error,1.0,-1.0);

					costs[iChannelOrder] = Blas_Dot_Prod(error,error);
					computedObservationsVector[iChannelOrder] = computedObservations;
				}

				int iLeastCost;
				Util::Min(costs,iLeastCost);

				particleCandidates[iCandidate].fromParticle = iParticle;
				particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;
                particleCandidates[iCandidate].iBestChannelOrder = iLeastCost;
				particleCandidates[iCandidate].weight = processedParticle->GetWeight()*StatUtil::NormalPdf(observations.col(iObservationToBeProcessed),computedObservationsVector[iLeastCost],noiseVariances[iObservationToBeProcessed]);

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
			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			// sampled symbols are copied into the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.col(_maxOrder-1));

			for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				// channel matrix is estimated by means of the particle channel estimator
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed,processedParticle->GetChannelMatrixEstimator(iChannelOrder)->NextMatrix(observations.col(iObservationToBeProcessed),particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix(rAllSymbolRows,tRange(_maxOrder-_candidateOrders[iChannelOrder],_maxOrder-1)),noiseVariances[iObservationToBeProcessed]));
			}

			processedParticle->SetWeight(particleCandidates[indexesSelectedCandidates[iParticle]].weight);

		} // for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)

		_particleFilter->NormalizeWeights();
	}

	delete[] particleCandidates;
}

int UPSPBasedSMCAlgorithm::BestChannelOrderIndex(int iBestParticle)
{
// 	ParticleWithChannelEstimation *bestParticle = dynamic_cast <ParticleWithChannelEstimation *>(_particleFilter->GetParticle(iBestParticle));

// 	int iMaxChannelOrderAPP = 0;
// 	double maxChannelOrderAPP = bestParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(iMaxChannelOrderAPP);
//
// 	for(uint i=1;i<_candidateOrders.size();i++)
// 		if(bestParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(i) > maxChannelOrderAPP)
// 		{
// 			maxChannelOrderAPP = bestParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(i);
// 			iMaxChannelOrderAPP = i;
// 		}
//
// 	return iMaxChannelOrderAPP;
// 	return 0;

    cout << "El orden bueno es " << _particlesBestChannelOrders[iBestParticle] << endl;

    return _particlesBestChannelOrders[iBestParticle];
}
