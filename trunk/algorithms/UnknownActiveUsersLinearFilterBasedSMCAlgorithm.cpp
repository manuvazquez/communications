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
#include "UnknownActiveUsersLinearFilterBasedSMCAlgorithm.h"

UnknownActiveUsersLinearFilterBasedSMCAlgorithm::UnknownActiveUsersLinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet, int L, int Nr, int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, LinearDetector *linearDetector, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances): SMCAlgorithm(name, alphabet, L, Nr, N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances),_linearDetector(linearDetector->clone())
{
    _randomParticlesInitilization = true; 
}


UnknownActiveUsersLinearFilterBasedSMCAlgorithm::~UnknownActiveUsersLinearFilterBasedSMCAlgorithm()
{
    delete _linearDetector;
}


void UnknownActiveUsersLinearFilterBasedSMCAlgorithm::initializeParticles()
{
    ChannelMatrixEstimator *channelMatrixEstimatorClone;
    tVector channelMean = Util::toVector(_channelMatrixMean,rowwise);
    tMatrix channelCovariance = LaGenMatDouble::from_diag(Util::toVector(_channelMatrixVariances,rowwise));

    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
    {
        channelMatrixEstimatorClone = _channelEstimator->clone();
        if(_randomParticlesInitilization)
            channelMatrixEstimatorClone->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::RandMatrix(channelMean,channelCovariance),rowwise,_nOutputs));
        
        _particleFilter->addParticle(new ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers(1.0/(double)_particleFilter->capacity(),_nInputs,_iLastSymbolVectorToBeDetected,channelMatrixEstimatorClone,_linearDetector->clone()));

        if(_preamble.cols()!=0)
            _particleFilter->getParticle(iParticle)->setSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);
    }
}

void UnknownActiveUsersLinearFilterBasedSMCAlgorithm::process(const tMatrix& observations, vector< double > noiseVariances)
{
    int iParticle,iSampledSymbol,iAlphabet,iSampled;
    tRange rAll;
    double proposal,s2q,sumProb,likelihoodsProd;
    tMatrix forWeightUpdateNeededSymbols(_nInputs,_channelOrder+_d);
    
    int iOutput;
    tMatrix channelMatrixSample;
    tVector sampledVector(_nInputs);
    tMatrix symbolProb(_nInputs,_alphabet.length());
    tVector predictedNoiselessObservation(_nOutputs);

    for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
        // noise covariance needs to be constructed
        tMatrix noiseCovariance = LaGenMatDouble::zeros(_nOutputs,_nOutputs);
        for(iOutput=0;iOutput<_nOutputs;iOutput++)
            noiseCovariance(iOutput,iOutput) = noiseVariances[iObservationToBeProcessed];

        for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
        {
            ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers *> (_particleFilter->getParticle(iParticle));            
            
            channelMatrixSample = (dynamic_cast<KalmanEstimator *> (_particleFilter->getParticle(iParticle)->getChannelMatrixEstimator(_estimatorIndex)))->sampleFromPredictive();            

            tVector softEstimations;

            // the sampled channel matrix is used to obtain soft estimations of the transmitted symbols
            softEstimations =  processedParticle->GetLinearDetector(_estimatorIndex)->detect(observations.col(iObservationToBeProcessed),channelMatrixSample,noiseCovariance);

            // the evaluated proposal function (necessary for computing the weights) is initialized
            proposal = 1.0;

            // sampling
            for(iSampledSymbol=0;iSampledSymbol<_nInputs;iSampledSymbol++)
            {
                s2q = processedParticle->GetLinearDetector(_estimatorIndex)->nthSymbolVariance(iSampledSymbol);

                sumProb = 0.0;

                // the probability for each posible symbol alphabet is computed
                for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                {
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
                sampledVector(iSampledSymbol) = _alphabet[iSampled];

                proposal *= symbolProb(iSampledSymbol,iSampled);
            }

            // sampled symbol vector is stored for the corresponding particle
            processedParticle->setSymbolVector(iObservationToBeProcessed,sampledVector);

            // predictedNoiselessObservation = matricesToStack[iSmoothing] * stackedSymbolVector
            Blas_Mat_Vec_Mult(channelMatrixSample,sampledVector,predictedNoiselessObservation);

            likelihoodsProd = StatUtil::NormalPdf(observations.col(iObservationToBeProcessed),predictedNoiselessObservation,noiseVariances[iObservationToBeProcessed]);

            // the weight is updated
            processedParticle->setWeight((likelihoodsProd/proposal)*processedParticle->getWeight());

            // and the estimation of the channel matrix
            processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,
                                                processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),sampledVector,noiseVariances[iObservationToBeProcessed]));

        } // for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)

        _particleFilter->normalizeWeights();

        // if it's not the last time instant
        if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
            _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

