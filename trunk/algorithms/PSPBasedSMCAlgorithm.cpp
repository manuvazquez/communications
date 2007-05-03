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

#include <defines.h>

#ifdef EXPORT_REAL_DATA
	extern MIMOChannel *realChannel;
	extern tMatrix *realSymbols;
	extern Noise *realNoise;
#endif

PSPBasedSMCAlgorithm::PSPBasedSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances, double ARcoefficient): SMCAlgorithm(name, alphabet, L, N, K, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances),_ARcoefficient(ARcoefficient)
{
}

void PSPBasedSMCAlgorithm::InitializeParticles()
{
	#ifdef DEBUG_PSPBASEDSMCALGORITHM
		cout << "Principio InitializeParticles" << endl;
	#endif

	// we begin with only one particle
	_particleFilter->AddParticle(new ParticleWithChannelEstimation(1.0,_N,_K+_d,_channelEstimator->Clone()));
	_particleFilter->GetParticle(0)->SetSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);

	#ifdef DEBUG_PSPBASEDSMCALGORITHM
		cout << "Final InitializeParticles" << endl;
	#endif
}

void PSPBasedSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	#ifdef DEBUG_PSPBASEDSMCALGORITHM
		cout << "Al principio de process" << endl;
	#endif

	uint nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);
	tRange rmMinus1FirstColumns(0,_m-2),rAllSymbolRows(0,_N-1);
	vector<tSymbol> testedVector(_N);
	tVector computedObservations(_L);

	typedef struct{
		int fromParticle;
		tMatrix symbolVectorsMatrix;
		double weight;
	}tParticleCandidate;

	tParticleCandidate *particleCandidates = new tParticleCandidate[_particleFilter->Nparticles()*nSymbolVectors] ;

	#ifdef DEBUG_PSPBASEDSMCALGORITHM
		cout << "Despues del new" << endl;
		cout << "_startDetectionTime es " << _startDetectionTime << endl;
	#endif

    // "symbolVectorsMatrix" will contain all the symbols involved in the current observation
    tMatrix symbolVectorsMatrix(_N,_m);
    tVector symbolsVector(_Nm);

    int lastSymbolVectorStart = _Nm - _N;

	// at first, there is only one particle
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K+_d;iObservationToBeProcessed++)
	{
		tRange rmMinus1PrecedentColumns(iObservationToBeProcessed-_m+1,iObservationToBeProcessed-1);

		// it keeps track of the place where a new tParticleCandidate will be stored within the vector
		int iCandidate = 0;

		#ifdef DEBUG_PSPBASEDSMCALGORITHM
			cout << "Antes de expandir las partículas los pesos son" << endl << _particleFilter->GetWeightsVector();
		#endif

		// the candidates from all the particles are generated
		for(int iParticle=0;iParticle<_particleFilter->NactualParticles();iParticle++)
		{
			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			symbolVectorsMatrix(rAllSymbolRows,rmMinus1FirstColumns).inject(processedParticle->GetSymbolVectors(rmMinus1PrecedentColumns));

			symbolsVector = Util::ToVector(symbolVectorsMatrix,columnwise);

			#ifdef DEBUG_PSPBASEDSMCALGORITHM
// 				cout << "Despues del inject" << endl;
			#endif

			tMatrix estimatedChannelMatrix = processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->LastEstimatedChannelMatrix();
			estimatedChannelMatrix *= _ARcoefficient;

			#ifdef DEBUG_PSPBASEDSMCALGORITHM
				cout << "La estimación de canal es" << endl << estimatedChannelMatrix;
				#ifdef EXPORT_REAL_DATA
					tMatrix verdaderoCanal = (*realChannel)[iObservationToBeProcessed];
					cout << "El verdadero canal es " << endl << verdaderoCanal;
					tMatrix verdaderosSimbolos = (*realSymbols)(tRange(0,_N-1),tRange(iObservationToBeProcessed-_m+1,iObservationToBeProcessed));
// 					cout << "Los verdaderos simbolos son" << endl << verdaderosSimbolos;
					tVector verdaderoRuido = (*realNoise)[iObservationToBeProcessed];
					tVector observacionCalculada(_L);
					Blas_Mat_Vec_Mult(verdaderoCanal,Util::ToVector(verdaderosSimbolos,columnwise),observacionCalculada);
					Util::Add(observacionCalculada,verdaderoRuido,observacionCalculada);
	// 				cout << "La observacion calculada" << endl << observacionCalculada;
				#endif
// 				cout << "La que tocaría" << endl << observations.col(iObservationToBeProcessed);
// 				cout << "Todas las observaciones" << endl << observations;
				#ifdef EXPORT_REAL_DATA
	// 				cout << "Todos los simbolos" << endl << (*realSymbols);
				#endif
				if(iParticle==1)
				{
					cout << "Los simbolos hasta ahora detectados" << endl << processedParticle->GetAllSymbolVectors()(tRange(0,_N-1),tRange(0,iObservationToBeProcessed-1));
					#ifdef EXPORT_REAL_DATA
						cout << "Los simbolos de verdad" << endl << (*realSymbols)(tRange(0,_N-1),tRange(0,iObservationToBeProcessed-1));
					#endif
				}
			#endif

			for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.IntToSymbolsArray(iTestedVector,testedVector);

// 				// current tested vector is copied in the m-th position
// 				for(int k=0;k<_N;k++)
// 					symbolVectorsMatrix(k,_m-1) = testedVector[k];

				// current tested vector is copied in the m-th position
				for(int k=0;k<_N;k++)
					symbolVectorsMatrix(k,_m-1) = symbolsVector(lastSymbolVectorStart+k) = testedVector[k];

// 				// computedObservations = estimatedChannelMatrix * symbolVectorsMatrix(:)
// 				Blas_Mat_Vec_Mult(estimatedChannelMatrix,Util::ToVector(symbolVectorsMatrix,columnwise),computedObservations);

				// computedObservations = estimatedChannelMatrix * symbolVectorsMatrix(:)
				Blas_Mat_Vec_Mult(estimatedChannelMatrix,symbolsVector,computedObservations);

				particleCandidates[iCandidate].fromParticle = iParticle;
				particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;
				particleCandidates[iCandidate].weight = processedParticle->GetWeight()*StatUtil::NormalPdf(observations.col(iObservationToBeProcessed),computedObservations,noiseVariances[iObservationToBeProcessed]);

				#ifdef DEBUG
					cout << "El peso actualizado para iTestedVector = " << iTestedVector << " es " << particleCandidates[iCandidate].weight << endl;
				#endif

				#ifdef DEBUG_PSPBASEDSMCALGORITHM
					if(iParticle==1)
						cout << "El peso actualizado es " << particleCandidates[iCandidate].weight << endl;
				#endif

				iCandidate++;
			} // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

			#ifdef DEBUG
				cout << "Una tecla..."; getchar();
			#endif

		} // for(int iParticle=0;iParticle<_particleFilter->NactualParticles();iParticle++)

		#ifdef DEBUG_PSPBASEDSMCALGORITHM
			cout << "Expandidas las particulas" << endl;
		#endif

		// a vector of size the number of generated candidates is declared...
		tVector weights(iCandidate);

		// ...to store their weights
		for(int i=0;i<iCandidate;i++)
			weights(i) = particleCandidates[i].weight;

		#ifdef DEBUG_PSPBASEDSMCALGORITHM
			cout << "Los pesos" << endl << weights;
// 			cout << "Una tecla..."; getchar();
		#endif

		// the candidates that are going to give particles are selected
		vector<int> indexesSelectedCandidates = _resamplingAlgorithm->ObtainIndexes(_particleFilter->Nparticles(),weights);

		// every survivor candidate is associated with an old particle
		vector<int> indexesParticles(indexesSelectedCandidates.size());
		for(uint i=0;i<indexesSelectedCandidates.size();i++)
			indexesParticles[i] = particleCandidates[indexesSelectedCandidates[i]].fromParticle;

		#ifdef DEBUG_PSPBASEDSMCALGORITHM
			cout << "Los indices de las particulas" << endl;
			Util::Print(indexesParticles);
		#endif

		// the chosen particles are kept without modification (yet)
		_particleFilter->KeepParticles(indexesParticles);

		#ifdef DEBUG_PSPBASEDSMCALGORITHM
			cout << "despues de keepParticles.." << endl;
		#endif

		// every surviving particle is modified according to what it says its corresponding candidate
		for(int iParticle=0;iParticle<_particleFilter->NactualParticles();iParticle++)
		{
			#ifdef DEBUG_PSPBASEDSMCALGORITHM
				cout << "iParticle " << iParticle << " iObservationToBeProcessed es " << iObservationToBeProcessed << endl;
				cout << particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix;
			#endif

			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			// sampled symbols are copied into the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.col(_m-1));

			#ifdef DEBUG_PSPBASEDSMCALGORITHM
// 				cout << "despues de SetSymbol" << endl;
			#endif

			// channel matrix is estimated by means of the particle channel estimator
			processedParticle->SetChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->NextMatrix(observations.col(iObservationToBeProcessed),particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix,noiseVariances[iObservationToBeProcessed]));

			#ifdef DEBUG_PSPBASEDSMCALGORITHM
// 				cout << "despues de SetChannel" << endl;
			#endif

			processedParticle->SetWeight(particleCandidates[indexesSelectedCandidates[iParticle]].weight);

			#ifdef DEBUG_PSPBASEDSMCALGORITHM
// 				cout << "despues de SetWeight" << endl;
			#endif
		} // for(int iParticle=0;iParticle<_particleFilter->NactualParticles();iParticle++)

		_particleFilter->NormalizeWeights();
	}

	delete[] particleCandidates;
}

