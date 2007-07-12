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
#include "LinearFilterBasedSMCAlgorithm.h"
#include <MMSEDetector.h>
#include <DecorrelatorDetector.h>

// #define SUBSTRACT_CONTRIBUTION_FROM_KNOWN_SYMBOLS

#define DEBUG
#define MUESTREO
#define DEBUG_SIN_RESTAR


#ifdef DEBUG_SIN_RESTAR
extern MIMOChannel *realChannel;
extern tMatrix *realSymbols;
extern Noise *realNoise;
#endif

LinearFilterBasedSMCAlgorithm::LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int N, int K,int m,  ChannelMatrixEstimator *channelEstimator,LinearDetector *linearDetector,tMatrix preamble, int backwardsSmoothingLag, int smoothingLag, int nParticles,ResamplingAlgorithm *resamplingAlgorithm,const tMatrix &channelMatrixMean, const tMatrix &channelMatrixVariances,double ARcoefficient,double samplingVariance,double ARprocessVariance): SMCAlgorithm(name, alphabet, L, N, K,m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances)
,_linearDetector(linearDetector->Clone()),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_c(backwardsSmoothingLag)
{
}

LinearFilterBasedSMCAlgorithm::LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int N, int K,int m,tMatrix preamble, int smoothingLag, ParticleFilter *particleFilter, ResamplingAlgorithm *resamplingAlgorithm,double ARcoefficient,double samplingVariance, double ARprocessVariance): SMCAlgorithm(name, alphabet, L, N, K,m, preamble, smoothingLag, particleFilter, resamplingAlgorithm)
,_linearDetector(NULL),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_c(0)
{
}

LinearFilterBasedSMCAlgorithm::~LinearFilterBasedSMCAlgorithm()
{
    delete _linearDetector;
}

void LinearFilterBasedSMCAlgorithm::InitializeParticles()
{
    tRange rPreamble(0,_preamble.cols()-1);

	// memory is reserved
	for(int iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)
	{
		_particleFilter->AddParticle(new ParticleWithChannelEstimationAndLinearDetection(1.0/(double)_particleFilter->Capacity(),_N,_K,_channelEstimator->Clone(),_linearDetector->Clone()));

        _particleFilter->GetParticle(iParticle)->SetSymbolVectors(rPreamble,_preamble);
	}
}

void LinearFilterBasedSMCAlgorithm::Process(const tMatrix &observations, vector< double > noiseVariances)
{
	int iParticle,iSmoothing,iRow,iSampledSymbol,iAlphabet,iSampled;
	vector<tMatrix> matricesToStack(_c+_d+1,tMatrix(_L,_Nm));
	tRange rAll,rAllSymbolRows(0,_N-1);
	tRange range0mMinus2(0,_m-2),rSampledSymbolVectors(_m-1,_m+_d-1);
	tRange rFirstmSymbolVectors(0,_m-1);
	tVector sampledVector(_N),sampledSmoothingVector(_N*(_d+1));
	double proposal,s2q,sumProb,likelihoodsProd;
	tMatrix s2qAux(_L*(_c+_d+1),_L*(_c+_d+1)),symbolProb(_N*(_d+1),_alphabet.Length());
	tVector s2qAuxFilter(_L*(_c+_d+1));
	tMatrix forWeightUpdateNeededSymbols(_N,_m+_d);
	tVector predictedNoiselessObservation(_L);

#ifdef MUESTREO
	tMatrix symbolProbAntes(_N*(_d+1),_alphabet.Length());
	double sumProbAntes;
#endif
#ifdef DEBUG_RESTANDO
	MMSEDetector mmseDetector(_L*(_c+_d+1),_N*(_d+1),_alphabet.Variance(),_N*(_d+1));
	DecorrelatorDetector decorrelatorDetector(_L*(_c+_d+1),_N*(_d+1),_alphabet.Variance());
#endif

#ifdef DEBUG_SIN_RESTAR
	MMSEDetector mmseDetector(_L*(_c+_d+1),_N*(_c+_m+_d),_alphabet.Variance(),_N*(_d+1));
	DecorrelatorDetector decorrelatorDetector(_L*(_c+_d+1),_N*(_d+1),_alphabet.Variance());
#endif

#ifdef SUBSTRACT_CONTRIBUTION_FROM_KNOWN_SYMBOLS
	tRange rColumnsStackedChannelMatrixMinus((_c+_m-1)*_N,(_c+_m+_d)*_N-1);
	if(_linearDetector->ChannelMatrixCols() != _N*(_d+1))
		throw RuntimeException("LinearFilterBasedSMCAlgorithm::Process: SUBSTRACT_CONTRIBUTION_FROM_KNOWN_SYMBOLS is defined. It shouldn't because it's not compatible with the current linear detector.");
#endif

	// already detected symbol vectors involved in the current detection
	tRange rmMinus1AlreadyDetectedSymbolVectors(_startDetectionTime-_m+1,_startDetectionTime-1);

	// observation matrix columns that are involved in the smoothing
	tRange rSmoothingRange(_startDetectionTime-_c,_startDetectionTime+_d);

#ifdef SUBSTRACT_CONTRIBUTION_FROM_KNOWN_SYMBOLS
	tRange rAlreadyDetectedSymbolVectors(_startDetectionTime-_c-_m+1,_startDetectionTime-1);
#endif

	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
#ifdef DEBUG
		cout << "----------------------------------------------------------------------" << endl;
#endif

		// the stacked observations vector
		tVector stackedObservations = Util::ToVector(observations(rAll,rSmoothingRange),columnwise);

		// stacked noise covariance needs to be constructed
		tMatrix stackedNoiseCovariance = LaGenMatDouble::zeros(_L*(_c+_d+1),_L*(_c+_d+1));
		for(iSmoothing=-_c;iSmoothing<=_d;iSmoothing++)
			for(iRow=0;iRow<_L;iRow++)
				stackedNoiseCovariance((iSmoothing+_c)*_L+iRow,(iSmoothing+_c)*_L+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];

		for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)
		{
			ParticleWithChannelEstimationAndLinearDetection *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetection *> (_particleFilter->GetParticle(iParticle));

			// already estimated channel matrices are stored in a vector in order to stack them
			for(iSmoothing=-_c;iSmoothing<0;iSmoothing++)
				matricesToStack[iSmoothing+_c] = processedParticle->GetChannelMatrix(_estimatorIndex,iObservationToBeProcessed+iSmoothing);

			// first of the predicted ones is obtained via a virtual method
			FillFirstEstimatedChannelMatrix(iParticle,matricesToStack[_c]);

			for(iSmoothing=_c+1;iSmoothing<=_c+_d;iSmoothing++)
				// matricesToStack[iSmoothing] = _ARcoefficient * matricesToStack[iSmoothing-1] + rand(_L,_Nm)*_ARprocessVariance
				Util::Add(matricesToStack[iSmoothing-1],StatUtil::RandnMatrix(_L,_Nm,0.0,_ARprocessVariance),matricesToStack[iSmoothing],_ARcoefficient,1.0);

			// matrices are stacked to give
			tMatrix stackedChannelMatrix = HsToStackedH(matricesToStack);

#ifdef DEBUG_SIN_RESTAR
			cout << "---------------------------------------------" << endl;
			vector<tMatrix> realMatricesToStack(_c+_d+1,tMatrix(_L,_Nm));
			for(iSmoothing=_c;iSmoothing<=_c+_d;iSmoothing++)
				realMatricesToStack[iSmoothing] = (*realChannel)[iObservationToBeProcessed+iSmoothing];

			tVector stackedNoise = Util::ToVector(realNoise->Range(iObservationToBeProcessed-_c,iObservationToBeProcessed+_d),columnwise);
			cout << "stackedObservations" << endl << stackedObservations;
			cout << "stackedChannelMatrix" << endl << stackedChannelMatrix;
			cout << "stackedNoise" << endl << stackedNoise;
			cout << "stackedNoiseCovariance" << endl << stackedNoiseCovariance;
			tMatrix realStackedChannelMatrix = HsToStackedH(realMatricesToStack);
			cout << "stackedChannelMatrix de verdad" << endl << realStackedChannelMatrix;
			cout << "los simbolos verdaderos" << endl << Util::ToVector((*realSymbols)(rAll,tRange(iObservationToBeProcessed,iObservationToBeProcessed+_d)),columnwise);
			cout << "soft estimations conociendo el canal" << endl << mmseDetector.Detect(stackedObservations,realStackedChannelMatrix,stackedNoiseCovariance);
#endif

#ifdef DEBUG_RESTANDO
			vector<tMatrix> realMatricesToStack(_c+_d+1,tMatrix(_L,_Nm));

			for(iSmoothing=_c;iSmoothing<=_c+_d;iSmoothing++)
				realMatricesToStack[iSmoothing] = (*realChannel)[iObservationToBeProcessed+iSmoothing];

			tMatrix realStackedChannelMatrix = HsToStackedH(realMatricesToStack);
			tVector realTransformedObservations = SubstractKnownSymbolsContribution(realMatricesToStack,_m,_c,_d,stackedObservations,(*realSymbols)(rAll,rAlreadyDetectedSymbolVectors));
			tVector stackedNoise = Util::ToVector(realNoise->Range(iObservationToBeProcessed-_c,iObservationToBeProcessed+_d),columnwise);
			tVector allInvolvedSymbols = Util::ToVector((*realSymbols)(rAll,tRange(iObservationToBeProcessed-_c-_m+1,iObservationToBeProcessed+_d)),columnwise);
			cout << "todos los símbolos involucrados" << endl << allInvolvedSymbols;
			cout << "ruido apilado" << endl << stackedNoise << endl;
			tVector observacionesCalculadas = stackedNoise;
			Blas_Mat_Vec_Mult(realStackedChannelMatrix,allInvolvedSymbols,observacionesCalculadas,1.0,1.0);
// 			cout << "las observaciones calculadas aquí dentro son" << endl << observacionesCalculadas;
// 			cout << "las que recibi" << endl << stackedObservations;
			realStackedChannelMatrix = realStackedChannelMatrix(rAll,rColumnsStackedChannelMatrixMinus);
			tVector unknownSymbols = Util::ToVector((*realSymbols)(rAll,tRange(iObservationToBeProcessed,iObservationToBeProcessed+_d)),columnwise);
			tVector realObservacionesTransformadas = stackedNoise;
			Blas_Mat_Vec_Mult(realStackedChannelMatrix,unknownSymbols,realObservacionesTransformadas,1.0,1.0);
// 			cout << "las observaciones transformadas de verdad" << endl << realObservacionesTransformadas;
			cout << "las observaciones transformadas de verdad" << endl << realTransformedObservations;
			cout << "la matriz de canal de verdad ya recortada" << endl << realStackedChannelMatrix;
			cout << "los simbolos desconocidos" << endl << unknownSymbols;
			cout << "la covarianza del ruido" << endl << stackedNoiseCovariance;
#endif

			// the estimated stacked channel matrix is used to obtain soft estimations of the transmitted symbols
#ifdef SUBSTRACT_CONTRIBUTION_FROM_KNOWN_SYMBOLS
			stackedChannelMatrix = stackedChannelMatrix(rAll,rColumnsStackedChannelMatrixMinus);
			tVector transformedObservations = SubstractKnownSymbolsContribution(matricesToStack,_m,_c,_d,stackedObservations,processedParticle->GetSymbolVectors(rAlreadyDetectedSymbolVectors));
#ifdef DEBUG_RESTANDO
			cout << "observaciones transformadas" << endl << transformedObservations;
#endif
			tVector softEstimations =  processedParticle->GetLinearDetector(_estimatorIndex)->Detect(transformedObservations,stackedChannelMatrix,stackedNoiseCovariance);
#else
			tVector softEstimations =  processedParticle->GetLinearDetector(_estimatorIndex)->Detect(stackedObservations,stackedChannelMatrix,stackedNoiseCovariance);
#endif

			tMatrix filter = processedParticle->GetLinearDetector(_estimatorIndex)->ComputedFilter();

#ifdef DEBUG_RESTANDO
			cout << "stackedChannelMatrix:" << endl << stackedChannelMatrix;
			cout << "filter:" << endl << filter;
// 			cout << "true symbols:" << endl << realSymbols->col(iObservationToBeProcessed);
			cout << "soft estimations:" << endl << softEstimations;
#endif

#ifdef DEBUG_SIN_RESTAR
			cout << "soft estimations:" << endl << softEstimations;
#endif

#ifdef DEBUG_RESTANDO
			cout << "Detect para la de verdad" << endl;
			cout << "soft estimations conociendo el canal" << endl << mmseDetector.Detect(realTransformedObservations,realStackedChannelMatrix,stackedNoiseCovariance);
#endif
            // the evaluated proposal function (necessary for computing the weights) is initialized
            proposal = 1.0;

			// sampling
			for(iSampledSymbol=0;iSampledSymbol<(_N*(_d+1));iSampledSymbol++)
			{
				s2q = processedParticle->GetLinearDetector(_estimatorIndex)->NthSymbolVariance(iSampledSymbol);

#ifdef MUESTREO
				cout << "******** s2q para el símbolo " << iSampledSymbol << ": " << s2q << " ganancia = " << processedParticle->GetLinearDetector(_estimatorIndex)->NthSymbolGain(iSampledSymbol) << " ************" << endl;
// 				cout << "s2q para el símbolo " << iSampledSymbol << ": " << s2q << endl;
#endif

				sumProb = 0.0;
#ifdef MUESTREO
				sumProbAntes = 0.0;
#endif
				// the probability for each posible symbol alphabet is computed
				for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
				{
// 					symbolProb(iSampledSymbol,iAlphabet) = StatUtil::NormalPdf(softEstimations(iSampledSymbol),_alphabet[iAlphabet],s2q);
					symbolProb(iSampledSymbol,iAlphabet) = StatUtil::NormalPdf(softEstimations(iSampledSymbol),processedParticle->GetLinearDetector(_estimatorIndex)->NthSymbolGain(iSampledSymbol)*_alphabet[iAlphabet],s2q);

#ifdef MUESTREO
					symbolProbAntes(iSampledSymbol,iAlphabet) = StatUtil::NormalPdf(softEstimations(iSampledSymbol),_alphabet[iAlphabet],s2q);
// 					cout << "la probabilidad calculada ahora: " << StatUtil::NormalPdf(softEstimations(iSampledSymbol),processedParticle->GetLinearDetector(_estimatorIndex)->NthSymbolGain(iSampledSymbol)*_alphabet[iAlphabet],s2q) << " antes: " << StatUtil::NormalPdf(softEstimations(iSampledSymbol),_alphabet[iAlphabet],s2q) << endl;
#endif

					// the computed pdf is accumulated for normalizing purposes
					sumProb += symbolProb(iSampledSymbol,iAlphabet);

#ifdef MUESTREO
					sumProbAntes += symbolProbAntes(iSampledSymbol,iAlphabet);
#endif
				}

				try {
					for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
					{
						symbolProb(iSampledSymbol,iAlphabet) /= sumProb;
#ifdef MUESTREO
						symbolProbAntes(iSampledSymbol,iAlphabet) /= sumProbAntes;
#endif
					}
				}catch(exception e){
					cout << "The sum of the probabilities is null." << endl;
					for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
						symbolProb(iSampledSymbol,iAlphabet) = 1.0/double(_alphabet.Length());
				}

#ifdef MUESTREO
				cout << "Ahora voy a muestrear de" << endl << symbolProb.row(iSampledSymbol);
// 				<< "antes de" << endl << symbolProbAntes.row(iSampledSymbol);
#endif

				iSampled = StatUtil::Discrete_rnd(symbolProb.row(iSampledSymbol));
				sampledSmoothingVector(iSampledSymbol) = _alphabet[iSampled];

				proposal *= symbolProb(iSampledSymbol,iSampled);
			}

#ifdef DEBUG_RESTANDO
			cout << "el vector de símbolos muestreado:" << endl << sampledSmoothingVector;
			if(sampledSmoothingVector(0)!=(*realSymbols)(0,iObservationToBeProcessed) || sampledSmoothingVector(1)!=(*realSymbols)(1,iObservationToBeProcessed))
			{
				cout << "Muestreó mal (iObservationToBeProcessed = " << iObservationToBeProcessed << ",iParticle = " << iParticle << ")" << endl;
				cout << "Una tecla..."; getchar();
			}
#endif

#ifdef DEBUG_SIN_RESTAR
			cout << "el vector de símbolos muestreado:" << endl << sampledSmoothingVector;
			if(sampledSmoothingVector(0)!=(*realSymbols)(0,iObservationToBeProcessed) || sampledSmoothingVector(1)!=(*realSymbols)(1,iObservationToBeProcessed))
			{
				cout << "Muestreó mal (iObservationToBeProcessed = " << iObservationToBeProcessed << ",iParticle = " << iParticle << ")" << endl;
				cout << "Una tecla..."; getchar();
			}
#endif

			// sampled symbol vector is stored for the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledSmoothingVector(rAllSymbolRows));

			// all the symbol vectors involved in the smoothing are kept in "forWeightUpdateNeededSymbols"
			// i) the already known:
			forWeightUpdateNeededSymbols(rAll,range0mMinus2).inject(processedParticle->GetSymbolVectors(rmMinus1AlreadyDetectedSymbolVectors));

			// ii) the just sampled
			forWeightUpdateNeededSymbols(rAll,rSampledSymbolVectors).inject(Util::ToMatrix(sampledSmoothingVector,columnwise,_N));

			likelihoodsProd = SmoothedLikelihood(matricesToStack,forWeightUpdateNeededSymbols,processedParticle,iObservationToBeProcessed,observations,noiseVariances);

			// the weight is updated
			processedParticle->SetWeight((likelihoodsProd/proposal)*processedParticle->GetWeight());

			// and the estimation of the channel matrix
			processedParticle->SetChannelMatrix(_estimatorIndex,iObservationToBeProcessed,
			                                    processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->NextMatrix(observations.col(iObservationToBeProcessed),
				                                    forWeightUpdateNeededSymbols(rAll,rFirstmSymbolVectors),noiseVariances[iObservationToBeProcessed]));

		} // for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)

		_particleFilter->NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
            _resamplingAlgorithm->ResampleWhenNecessary(_particleFilter);

		rmMinus1AlreadyDetectedSymbolVectors = rmMinus1AlreadyDetectedSymbolVectors + 1;
		rSmoothingRange = rSmoothingRange + 1;
#ifdef SUBSTRACT_CONTRIBUTION_FROM_KNOWN_SYMBOLS
		rAlreadyDetectedSymbolVectors = rAlreadyDetectedSymbolVectors + 1;
#endif

	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
}

vector<tMatrix> LinearFilterBasedSMCAlgorithm::ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
	int lengthSequenceToProcess = _preamble.cols() + trainingSequence.cols();
	tRange allObservationRows(0,_L-1);

	for(int i=_preamble.cols();i<lengthSequenceToProcess;i++)
	{
		tRange rSmoothingRange(i-_c,i+_d);
		tVector stackedObservationsVector = Util::ToVector(observations(allObservationRows,rSmoothingRange),columnwise);
		_linearDetector->StateStep(stackedObservationsVector);
	}

	return SMCAlgorithm::ProcessTrainingSequence(observations,noiseVariances,trainingSequence);
}
