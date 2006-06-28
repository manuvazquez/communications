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

LinearFilterBasedSMCAlgorithm::LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int N, int K,  ChannelMatrixEstimator *channelEstimator,LinearDetector *linearDetector,tMatrix preamble, int smoothingLag, int nParticles, ResamplingCriterion resamplingCriterion, StdResamplingAlgorithm resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance): SMCAlgorithm(name, alphabet, L, N, K,  channelEstimator, preamble, smoothingLag, nParticles, resamplingCriterion, resamplingAlgorithm)
,_linearDetector(linearDetector->Clone()),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance)
{
}


LinearFilterBasedSMCAlgorithm::~LinearFilterBasedSMCAlgorithm()
{
    delete _linearDetector;
}

void LinearFilterBasedSMCAlgorithm::InitializeParticles()
{
	// memory is reserved
	for(int iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
	{
		_particleFilter.SetParticle(new ParticleWithChannelEstimationAndLinearDetection(1.0/(double)_particleFilter.Nparticles(),_N,_K,_channelEstimator->Clone(),_linearDetector->Clone()),iParticle);
	}
}

void LinearFilterBasedSMCAlgorithm::Process(const tMatrix &observations, vector< double > noiseVariances)
{
	int iParticle,iSmoothing,iRow,iSampledSymbol,iAlphabet,iSampled;
	vector<tMatrix> matricesToStack(_d+1,tMatrix(_L,_Nm));
	tRange allObservationsRows(0,_L-1),allSymbolRows(0,_N-1);
	tRange range0mMinus2(0,_m-2),rSampledSymbolVectors(_m-1,_m+_d-1);
	tRange rFirstmSymbolVectors(0,_m-1);
	tVector sampledVector(_N),sampledSmoothingVector(_N*(_d+1));
	double proposal,s2q,sumProb,likelihoodsProd;
	tMatrix s2qAux(_L*(_d+1),_L*(_d+1)),symbolProb(_N*(_d+1),_alphabet.Length());
	tVector s2qAuxFilter(_L*(_d+1));
	tMatrix forWeightUpdateNeededSymbols(_N,_m+_d);
	tMatrix noiseCovariances[_d+1];
	tVector predictedNoiselessObservation(_L);

	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
		// already detected symbol vectors involved in the current detection
		tRange alreadyDetectedSymbolVectors(iObservationToBeProcessed-_m+1,iObservationToBeProcessed-1);

		// observation matrix columns that are involved in the smoothing
		tRange smoothingRange(iObservationToBeProcessed,iObservationToBeProcessed+_d);

		// the stacked observations vector
		tVector stackedObservations = Util::ToVector(observations(allObservationsRows,smoothingRange),columnwise);

		// stacked noise covariance needs to be constructed
		tMatrix stackedNoiseCovariance = LaGenMatDouble::zeros(_L*(_d+1),_L*(_d+1));
		for(iSmoothing=0;iSmoothing<_d+1;iSmoothing++)
			for(iRow=0;iRow<_L;iRow++)
				stackedNoiseCovariance(iSmoothing*_L+iRow,iSmoothing*_L+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];

		// required noise covariances are computed from the noise variances
		for(iSmoothing=0;iSmoothing<=_d;iSmoothing++)
		{
			noiseCovariances[iSmoothing] = LaGenMatDouble::eye(_L);
			noiseCovariances[iSmoothing] *= noiseVariances[iObservationToBeProcessed+iSmoothing];
		}

		for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
		{
			ParticleWithChannelEstimation *processedParticle = _particleFilter.GetParticle(iParticle);
		
// 			// predicted channel matrices are stored in a vector in order to stack them

			// matricesToStack[0] = _ARcoefficient * <lastEstimatedChannelMatrix> + randn(_L,_Nm)*_samplingVariance
			Util::Add((processedParticle->GetChannelMatrixEstimator())->LastEstimatedChannelMatrix(),StatUtil::RandnMatrix(_L,_Nm,0.0,_samplingVariance),matricesToStack[0],_ARcoefficient,1.0);
// 			Util::Add((_particlesChannelMatrixEstimators[iParticle])->LastEstimatedChannelMatrix(),StatUtil::RandnMatrix(_L,_Nm,0.0,_samplingVariance),matricesToStack[0],_ARcoefficient,1.0);

			for(iSmoothing=1;iSmoothing<=_d;iSmoothing++)
			{
				// matricesToStack[iSmoothing] = _ARcoefficient * matricesToStack[iSmoothing-1] + rand(_L,_Nm)*_ARprocessVariance
				Util::Add(matricesToStack[iSmoothing-1],StatUtil::RandnMatrix(_L,_Nm,0.0,_ARprocessVariance),matricesToStack[iSmoothing],_ARcoefficient,1.0);
			}

			// matrices are stacked to give
			tMatrix stackedChannelMatrix = HsToStackedH(matricesToStack);

			// the estimated stacked channel matrix is used to obtain soft estimations
			// of the transmitted symbols
			tVector softEstimations =  (dynamic_cast <ParticleWithChannelEstimationAndLinearDetection *> (processedParticle)->GetLinearDetector())->Detect(stackedObservations,stackedChannelMatrix);
// 			tVector softEstimations =  (_particlesLinearDetectors[iParticle])->Detect(stackedObservations,stackedChannelMatrix);

// 			tMatrix filter = _particlesLinearDetectors[iParticle]->ComputedFilter();
			tMatrix filter = (dynamic_cast <ParticleWithChannelEstimationAndLinearDetection *> (processedParticle)->GetLinearDetector())->ComputedFilter();

			// operations needed to computed the sampling variance

			//s2qAux = _alphabet.Variance() * stackedChannelMatrix * stackedChannelMatrix^H
			Blas_Mat_Mat_Trans_Mult(stackedChannelMatrix,stackedChannelMatrix,s2qAux,_alphabet.Variance());

			// s2qAux = s2qAux + stackedNoiseCovariance
			Util::Add(s2qAux,stackedNoiseCovariance,s2qAux);

			proposal = 1.0;
			
			// sampling
			for(iSampledSymbol=0;iSampledSymbol<(_N*(_d+1));iSampledSymbol++)
			{
				// s2qAuxFilter = s2qAux * filter.col(iSampledSymbol)
				Blas_Mat_Vec_Mult(s2qAux,filter.col(iSampledSymbol),s2qAuxFilter);


				s2q = _alphabet.Variance()*(1.0 - 2.0*Blas_Dot_Prod(filter.col(iSampledSymbol),stackedChannelMatrix.col(iSampledSymbol))) + Blas_Dot_Prod(filter.col(iSampledSymbol),s2qAuxFilter);
// 				cout << "La varianza calculada: " << s2q << " t = " << iObservationToBeProcessed << endl;

				double sumProb = 0.0;
				// the probability for each posible symbol alphabet is computed
				for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
				{
					symbolProb(iSampledSymbol,iAlphabet) = StatUtil::NormalPdf(softEstimations(iSampledSymbol),_alphabet[iAlphabet],s2q);

					// the computed pdf is accumulated for normalizing purposes
					sumProb += symbolProb(iSampledSymbol,iAlphabet);
				}

				try {
					for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
						symbolProb(iSampledSymbol,iAlphabet) /= sumProb;
				}catch(exception){
					cout << "The sum of the probabilities is null." << endl;
					for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
						symbolProb(iSampledSymbol,iAlphabet) = 0.5;
				}
				
				int iSampled = (StatUtil::Discrete_rnd(1,symbolProb.row(iSampledSymbol)))[0];
				sampledSmoothingVector(iSampledSymbol) = _alphabet[iSampled];

				proposal *= symbolProb(iSampledSymbol,iSampled);
			}

			// sampled symbol vector is stored for the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledSmoothingVector(allSymbolRows));

			// all the symbol vectors involved in the smoothing are kept in "forWeightUpdateNeededSymbols" 
			// i) the already known:
			forWeightUpdateNeededSymbols(allSymbolRows,range0mMinus2).inject(processedParticle->GetSymbolVectors(alreadyDetectedSymbolVectors));
			
			// ii) the just sampled
			forWeightUpdateNeededSymbols(allSymbolRows,rSampledSymbolVectors).inject(Util::ToMatrix(sampledSmoothingVector,columnwise,_N));

			likelihoodsProd = 1.0;

			for(iSmoothing=0;iSmoothing<=_d;iSmoothing++)
			{
				tRange rSymbolVectors(iSmoothing,iSmoothing+_m-1);
				tVector stackedSymbolVector = Util::ToVector(forWeightUpdateNeededSymbols(allSymbolRows,rSymbolVectors),columnwise);

				// predictedNoiselessObservation = matricesToStack[iSmoothing] * stackedSymbolVector
				Blas_Mat_Vec_Mult(matricesToStack[iSmoothing],stackedSymbolVector,predictedNoiselessObservation);

				likelihoodsProd *= StatUtil::NormalPdf(observations.col(iObservationToBeProcessed+iSmoothing),predictedNoiselessObservation,noiseCovariances[iSmoothing]);
			}

			// the weight is updated
			processedParticle->SetWeight((likelihoodsProd/proposal)*processedParticle->GetWeight());

			// and the estimation of the channel matrix
			processedParticle->SetChannelMatrix(iObservationToBeProcessed, (processedParticle->GetChannelMatrixEstimator())->NextMatrix(observations.col(iObservationToBeProcessed),forWeightUpdateNeededSymbols(allSymbolRows,rFirstmSymbolVectors),noiseVariances[iObservationToBeProcessed]));

		}
		_particleFilter.NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
			_particleFilter.Resampling();
	}
}

vector<tMatrix> LinearFilterBasedSMCAlgorithm::ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
	int lengthSequenceToProcess = _preamble.cols() + trainingSequence.cols();
	tRange allObservationRows(0,_L-1);

	for(int i=_m-1;i<lengthSequenceToProcess;i++)
	{
		tRange smoothingRange(i,i+_d);
		tVector stackedObservationsVector = Util::ToVector(observations(allObservationRows,smoothingRange),columnwise);
		_linearDetector->StateStep(stackedObservationsVector);
	}

// 	// the resultant linear detector is copied into each particle
// 	for(int iParticle=0;iParticle<_nParticles;iParticle++)
// 		_particlesLinearDetectors[iParticle] = _linearDetector->Clone();

	return SMCAlgorithm::ProcessTrainingSequence(observations,noiseVariances,trainingSequence);
}

