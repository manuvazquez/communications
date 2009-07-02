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
#include "PSPBasedSMCAlgorithm.h"

// #define DEBUG

#include <defines.h>

#ifdef EXPORT_REAL_DATA
	extern MIMOChannel *realChannel;
	extern tMatrix *realSymbols;
	extern Noise *realNoise;
#endif

PSPBasedSMCAlgorithm::PSPBasedSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int frameLength, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances, double ARcoefficient): SMCAlgorithm(name, alphabet, L, N, frameLength, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances),_ARcoefficient(ARcoefficient)
{
}

void PSPBasedSMCAlgorithm::InitializeParticles()
{
	// we begin with only one particle
	_particleFilter->AddParticle(new ParticleWithChannelEstimation(1.0,_nInputs,_K+_d,_channelEstimator->Clone()));
	_particleFilter->GetParticle(0)->SetSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);
}

void PSPBasedSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	uint nSymbolVectors = (int) pow((double)_alphabet.length(),(double)_nInputs);
	tRange rmMinus1FirstColumns(0,_channelOrder-2),rAllSymbolRows(0,_nInputs-1);
	vector<tSymbol> testedVector(_nInputs);
	tVector computedObservations(_nOutputs);
	double normConst;

	typedef struct{
		int fromParticle;
		tMatrix symbolVectorsMatrix;
		double weight;
	}tParticleCandidate;

	tParticleCandidate *particleCandidates = new tParticleCandidate[_particleFilter->Capacity()*nSymbolVectors] ;

    // "symbolVectorsMatrix" will contain all the symbols involved in the current observation
    tMatrix symbolVectorsMatrix(_nInputs,_channelOrder);
    tVector symbolsVector(_nInputsXchannelOrder);

    int lastSymbolVectorStart = _nInputsXchannelOrder - _nInputs;

	// at first, there is only one particle
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K+_d;iObservationToBeProcessed++)
	{
		tRange rmMinus1PrecedentColumns(iObservationToBeProcessed-_channelOrder+1,iObservationToBeProcessed-1);

		// it keeps track of the place where a new tParticleCandidate will be stored within the vector
		int iCandidate = 0;

		normConst = 0.0;

		// the candidates from all the particles are generated
		for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
		{
			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			symbolVectorsMatrix(rAllSymbolRows,rmMinus1FirstColumns).inject(processedParticle->GetSymbolVectors(rmMinus1PrecedentColumns));

			symbolsVector = Util::ToVector(symbolVectorsMatrix,columnwise);

			tMatrix estimatedChannelMatrix = processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->lastEstimatedChannelMatrix();
			estimatedChannelMatrix *= _ARcoefficient;

			for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.int2symbolsArray(iTestedVector,testedVector);

				// current tested vector is copied in the m-th position
				for(int k=0;k<_nInputs;k++)
					symbolVectorsMatrix(k,_channelOrder-1) = symbolsVector(lastSymbolVectorStart+k) = testedVector[k];

				// computedObservations = estimatedChannelMatrix * symbolVectorsMatrix(:)
				Blas_Mat_Vec_Mult(estimatedChannelMatrix,symbolsVector,computedObservations);

				particleCandidates[iCandidate].fromParticle = iParticle;
				particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;
				particleCandidates[iCandidate].weight = processedParticle->GetWeight()*StatUtil::NormalPdf(observations.col(iObservationToBeProcessed),computedObservations,noiseVariances[iObservationToBeProcessed]);
				normConst += particleCandidates[iCandidate].weight;

				iCandidate++;
			} // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

		} // for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)

		// a vector of size the number of generated candidates is declared...
		tVector weights(iCandidate);

		// ...to store their weights
		for(int i=0;i<iCandidate;i++)
			weights(i) = particleCandidates[i].weight/normConst;

#ifdef DEBUG
		cout << "La suma es " << Util::Sum(weights) << endl;
		cout << "Antes de llamar a _resamplingAlgorithm->ObtainIndexes" << endl;
#endif

		// the candidates that are going to give particles are selected
		vector<int> indexesSelectedCandidates = _resamplingAlgorithm->ObtainIndexes(_particleFilter->Capacity(),weights);

#ifdef DEBUG
		cout << "Despues de llamar a _resamplingAlgorithm->ObtainIndexes" << endl;
#endif

		// every survivor candidate is associated with an old particle
		vector<int> indexesParticles(indexesSelectedCandidates.size());
		for(uint i=0;i<indexesSelectedCandidates.size();i++)
			indexesParticles[i] = particleCandidates[indexesSelectedCandidates[i]].fromParticle;

		// the chosen particles are kept without modification (yet)
		_particleFilter->KeepParticles(indexesParticles);

		// every surviving particle is modified according to what it says its corresponding candidate
		for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
		{
			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			// sampled symbols are copied into the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.col(_channelOrder-1));

			// channel matrix is estimated by means of the particle channel estimator
			processedParticle->SetChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix,noiseVariances[iObservationToBeProcessed]));

			processedParticle->SetWeight(particleCandidates[indexesSelectedCandidates[iParticle]].weight);
		} // for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)

		_particleFilter->NormalizeWeights();
	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K+_d;iObservationToBeProcessed++)

	delete[] particleCandidates;
}

