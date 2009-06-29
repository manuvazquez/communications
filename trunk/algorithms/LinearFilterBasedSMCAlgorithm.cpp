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

LinearFilterBasedSMCAlgorithm::LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int N, int frameLength,int m,  ChannelMatrixEstimator *channelEstimator,LinearDetector *linearDetector,tMatrix preamble, int backwardsSmoothingLag, int SMCsmoothingLag, int forwardSmoothingLag, int nParticles,ResamplingAlgorithm *resamplingAlgorithm,const tMatrix &channelMatrixMean, const tMatrix &channelMatrixVariances,double ARcoefficient,double samplingVariance,double ARprocessVariance, bool substractContributionFromKnownSymbols): SMCAlgorithm(name, alphabet, L, N, frameLength,m, channelEstimator, preamble, SMCsmoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances)
,_linearDetector(linearDetector->Clone()),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_c(backwardsSmoothingLag),_e(forwardSmoothingLag),_substractContributionFromKnownSymbols(substractContributionFromKnownSymbols)
{
}

LinearFilterBasedSMCAlgorithm::LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int N, int frameLength,int m,tMatrix preamble, int SMCsmoothingLag, ParticleFilter *particleFilter, ResamplingAlgorithm *resamplingAlgorithm,double ARcoefficient,double samplingVariance, double ARprocessVariance, bool substractContributionFromKnownSymbols): SMCAlgorithm(name, alphabet, L, N, frameLength,m, preamble, SMCsmoothingLag, particleFilter, resamplingAlgorithm)
,_linearDetector(NULL),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_c(0),_e(SMCsmoothingLag),_substractContributionFromKnownSymbols(substractContributionFromKnownSymbols)
{
}

LinearFilterBasedSMCAlgorithm::~LinearFilterBasedSMCAlgorithm()
{
    delete _linearDetector;
}

void LinearFilterBasedSMCAlgorithm::InitializeParticles()
{
    tRange rPreamble(0,_preamble.cols()-1);

    ChannelMatrixEstimator *channelMatrixEstimatorClone;
    tVector channelMean = Util::ToVector(_channelMatrixMean,rowwise);
    tMatrix channelCovariance = LaGenMatDouble::from_diag(Util::ToVector(_channelMatrixVariances,rowwise));

    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)
    {
        channelMatrixEstimatorClone = _channelEstimator->Clone();
        if(_randomParticlesInitilization)
            channelMatrixEstimatorClone->SetFirstEstimatedChannelMatrix(Util::ToMatrix(StatUtil::RandMatrix(channelMean,channelCovariance),rowwise,_L));
        _particleFilter->AddParticle(new ParticleWithChannelEstimationAndLinearDetection(1.0/(double)_particleFilter->Capacity(),_N,_K,channelMatrixEstimatorClone,_linearDetector->Clone()));

        _particleFilter->GetParticle(iParticle)->SetSymbolVectors(rPreamble,_preamble);
    }
}

void LinearFilterBasedSMCAlgorithm::Process(const tMatrix &observations, vector< double > noiseVariances)
{
	int iParticle,iSmoothing,iRow,iSampledSymbol,iAlphabet,iSampled;
	vector<tMatrix> matricesToStack(_c+_e+1,tMatrix(_L,_Nm));
	tRange rAll,rNfirst(0,_N-1);
	tRange range0mMinus2(0,_m-2),rSampledSymbolVectors(_m-1,_m+_d-1);
	tRange rFirstmSymbolVectors(0,_m-1);
	tVector sampledVector(_N),sampledSmoothingVector(_N*(_d+1));
	double proposal,s2q,sumProb,likelihoodsProd;
	tMatrix s2qAux(_L*(_c+_d+1),_L*(_c+_d+1)),symbolProb(_N*(_d+1),_alphabet.length());
	tVector s2qAuxFilter(_L*(_c+_d+1));
	tMatrix forWeightUpdateNeededSymbols(_N,_m+_d);
	tVector predictedNoiselessObservation(_L);

    // just needed when substracting the contribution from the known symbols
    tRange rAlreadyDetectedSymbolVectors;

    if(_substractContributionFromKnownSymbols)
	   if(_linearDetector->ChannelMatrixCols() != _N*(_e+1))
		  throw RuntimeException("LinearFilterBasedSMCAlgorithm::Process: the algorithm is supposed to operate substracting the contribution of the known symbols but this is not compatible with the current linear detector.");

	// already detected symbol vectors involved in the current detection
	tRange rmMinus1AlreadyDetectedSymbolVectors(_startDetectionTime-_m+1,_startDetectionTime-1);

	// observation matrix columns that are involved in the smoothing
	tRange rSmoothingRange(_startDetectionTime-_c,_startDetectionTime+_e);

    if(_substractContributionFromKnownSymbols)
        rAlreadyDetectedSymbolVectors = tRange(_startDetectionTime-_c-_m+1,_startDetectionTime-1);

	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
		// the stacked observations vector
		tVector stackedObservations = Util::ToVector(observations(rAll,rSmoothingRange),columnwise);

		// stacked noise covariance needs to be constructed
		tMatrix stackedNoiseCovariance = LaGenMatDouble::zeros(_L*(_c+_e+1),_L*(_c+_e+1));
		for(iSmoothing=-_c;iSmoothing<=_e;iSmoothing++)
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

			for(iSmoothing=_c+1;iSmoothing<=_c+_e;iSmoothing++)
				// matricesToStack[iSmoothing] = _ARcoefficient * matricesToStack[iSmoothing-1] + rand(_L,_Nm)*_ARprocessVariance
				Util::Add(matricesToStack[iSmoothing-1],StatUtil::RandnMatrix(_L,_Nm,0.0,_ARprocessVariance),matricesToStack[iSmoothing],_ARcoefficient,1.0);

			// matrices are stacked to give
			tMatrix stackedChannelMatrix = HsToStackedH(matricesToStack);


            tVector softEstimations;

			// the estimated stacked channel matrix is used to obtain soft estimations of the transmitted symbols
            if(_substractContributionFromKnownSymbols)
            {
                softEstimations =  processedParticle->GetLinearDetector(_estimatorIndex)->Detect(
                        // transformed observations
                        SubstractKnownSymbolsContribution(matricesToStack,_m,_c,_e,stackedObservations,processedParticle->GetSymbolVectors(rAlreadyDetectedSymbolVectors)),
                        // only a part of the channel matrix is needed. The first range chooses all the stacked observation rows
                        stackedChannelMatrix(rAll,tRange((_c+_m-1)*_N,(_c+_m+_e)*_N-1)),
                        stackedNoiseCovariance);
            } else
            {
                softEstimations =  processedParticle->GetLinearDetector(_estimatorIndex)->Detect(stackedObservations,stackedChannelMatrix,stackedNoiseCovariance);
            }

            // the evaluated proposal function (necessary for computing the weights) is initialized
            proposal = 1.0;

			// sampling
			for(iSampledSymbol=0;iSampledSymbol<(_N*(_d+1));iSampledSymbol++)
			{
				s2q = processedParticle->GetLinearDetector(_estimatorIndex)->nthSymbolVariance(iSampledSymbol);

				sumProb = 0.0;

				// the probability for each posible symbol alphabet is computed
				for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
				{
// 					symbolProb(iSampledSymbol,iAlphabet) = StatUtil::NormalPdf(softEstimations(iSampledSymbol),_alphabet[iAlphabet],s2q);
					symbolProb(iSampledSymbol,iAlphabet) = StatUtil::NormalPdf(softEstimations(iSampledSymbol),processedParticle->GetLinearDetector(_estimatorIndex)->nthSymbolGain(iSampledSymbol)*_alphabet[iAlphabet],s2q);

					// the computed pdf is accumulated for normalizing purposes
					sumProb += symbolProb(iSampledSymbol,iAlphabet);
				}

				try {
					for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
					{
						symbolProb(iSampledSymbol,iAlphabet) /= sumProb;
					}
				}catch(exception e){
					cout << "The sum of the probabilities is null." << endl;
					for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
						symbolProb(iSampledSymbol,iAlphabet) = 1.0/double(_alphabet.length());
				}

				iSampled = StatUtil::discrete_rnd(symbolProb.row(iSampledSymbol));
				sampledSmoothingVector(iSampledSymbol) = _alphabet[iSampled];

				proposal *= symbolProb(iSampledSymbol,iSampled);
			}

			// sampled symbol vector is stored for the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledSmoothingVector(rNfirst));

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

        if(_substractContributionFromKnownSymbols)
		  rAlreadyDetectedSymbolVectors = rAlreadyDetectedSymbolVectors + 1;

	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
}

void LinearFilterBasedSMCAlgorithm::BeforeInitializingParticles(const tMatrix &observations, const tMatrix &trainingSequence)
{
    _linearDetector->StateStepsFromObservationsSequence(observations,_d,_preamble.cols(),_preamble.cols()+trainingSequence.cols());
}
