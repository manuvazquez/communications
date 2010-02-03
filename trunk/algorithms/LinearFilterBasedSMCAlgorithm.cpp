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

LinearFilterBasedSMCAlgorithm::LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m,  ChannelMatrixEstimator *channelEstimator,LinearDetector *linearDetector,MatrixXd preamble, int backwardsSmoothingLag, int SMCsmoothingLag, int forwardSmoothingLag, int nParticles,ResamplingAlgorithm *resamplingAlgorithm,const MatrixXd &channelMatrixMean, const MatrixXd &channelMatrixVariances,double ARcoefficient,double samplingVariance,double ARprocessVariance, bool substractContributionFromKnownSymbols): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected,m, channelEstimator, preamble, SMCsmoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances)
,_linearDetector(linearDetector->clone()),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_c(backwardsSmoothingLag),_e(forwardSmoothingLag),_substractContributionFromKnownSymbols(substractContributionFromKnownSymbols)
{
}

LinearFilterBasedSMCAlgorithm::LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m,MatrixXd preamble, int SMCsmoothingLag, ParticleFilter *particleFilter, ResamplingAlgorithm *resamplingAlgorithm,double ARcoefficient,double samplingVariance, double ARprocessVariance, bool substractContributionFromKnownSymbols): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected,m, preamble, SMCsmoothingLag, particleFilter, resamplingAlgorithm)
,_linearDetector(NULL),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_c(0),_e(SMCsmoothingLag),_substractContributionFromKnownSymbols(substractContributionFromKnownSymbols)
{
}

LinearFilterBasedSMCAlgorithm::~LinearFilterBasedSMCAlgorithm()
{
    delete _linearDetector;
}

void LinearFilterBasedSMCAlgorithm::initializeParticles()
{
    ChannelMatrixEstimator *channelMatrixEstimatorClone;
    
    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
    {
        channelMatrixEstimatorClone = _channelEstimator->clone();
        if(_randomParticlesInitilization)
            channelMatrixEstimatorClone->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::randnMatrix(_channelMean,_channelCovariance),rowwise,_nOutputs));
        _particleFilter->addParticle(new ParticleWithChannelEstimationAndLinearDetection(1.0/(double)_particleFilter->capacity(),_nInputs,_iLastSymbolVectorToBeDetected,channelMatrixEstimatorClone,_linearDetector->clone()));

        _particleFilter->getParticle(iParticle)->setSymbolVectors(0,_preamble.cols(),_preamble);
    }
}

// eigen
void LinearFilterBasedSMCAlgorithm::process(const MatrixXd &observations, vector<double> noiseVariances)
{
    int iParticle,iSmoothing,iRow,iSampledSymbol,iAlphabet,iSampled;
    vector<MatrixXd> matricesToStack(_c+_e+1,MatrixXd(_nOutputs,_nInputsXchannelOrder));
    VectorXd sampledVector(_nInputs),sampledSmoothingVector(_nInputs*(_d+1));
    double proposal,s2q,sumProb,likelihoodsProd;
    MatrixXd s2qAux(_nOutputs*(_c+_d+1),_nOutputs*(_c+_d+1)),symbolProb(_nInputs*(_d+1),_alphabet.length());
    VectorXd s2qAuxFilter(_nOutputs*(_c+_d+1));
    MatrixXd forWeightUpdateNeededSymbols(_nInputs,_channelOrder+_d);
    VectorXd predictedNoiselessObservation(_nOutputs);

    if(_substractContributionFromKnownSymbols)
       if(_linearDetector->channelMatrixcols() != _nInputs*(_e+1))
          throw RuntimeException("LinearFilterBasedSMCAlgorithm::process: the algorithm is supposed to operate substracting the contribution of the known symbols but this is not compatible with the current linear detector.");

    for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
        // the stacked observations vector
        VectorXd stackedObservations = Util::toVector(observations.block(0,iObservationToBeProcessed-_c,_nOutputs,_c+_e+1),columnwise);

        // stacked noise covariance needs to be constructed
        MatrixXd stackedNoiseCovariance = MatrixXd::Zero(_nOutputs*(_c+_e+1),_nOutputs*(_c+_e+1));
        for(iSmoothing=-_c;iSmoothing<=_e;iSmoothing++)
            for(iRow=0;iRow<_nOutputs;iRow++)
                stackedNoiseCovariance((iSmoothing+_c)*_nOutputs+iRow,(iSmoothing+_c)*_nOutputs+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];

        for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
        {
            // already estimated channel matrices are stored in a vector in order to stack them
            for(iSmoothing=-_c;iSmoothing<0;iSmoothing++)
                matricesToStack[iSmoothing+_c] = dynamic_cast <WithChannelEstimationParticleAddon *> (_particleFilter->getParticle(iParticle))->getChannelMatrix_eigen(_estimatorIndex,iObservationToBeProcessed+iSmoothing);

            // first of the predicted ones is obtained via a virtual method
            fillFirstEstimatedChannelMatrix(iParticle,matricesToStack[_c]);

            for(iSmoothing=_c+1;iSmoothing<=_c+_e;iSmoothing++)
                // matricesToStack[iSmoothing] = _ARcoefficient * matricesToStack[iSmoothing-1] + rand(_nOutputs,_nInputsXchannelOrder)*_ARprocessVariance
                matricesToStack[iSmoothing] = _ARcoefficient*matricesToStack[iSmoothing-1] + StatUtil::randnMatrix(_nOutputs,_nInputsXchannelOrder,0.0,_ARprocessVariance);

            // matrices are stacked to give
            MatrixXd stackedChannelMatrix = channelMatrices2stackedChannelMatrix(matricesToStack);

            VectorXd softEstimations;

            // the estimated stacked channel matrix is used to obtain soft estimations of the transmitted symbols
            if(_substractContributionFromKnownSymbols)
            {
                // transformed observations
                softEstimations =  dynamic_cast <WithLinearDetectionParticleAddon *> (_particleFilter->getParticle(iParticle))->getLinearDetector(_estimatorIndex)->detect(substractKnownSymbolsContribution(matricesToStack,_channelOrder,_c,_e,stackedObservations,_particleFilter->getParticle(iParticle)->getSymbolVectors(iObservationToBeProcessed-_c-_channelOrder+1,iObservationToBeProcessed-1)),stackedChannelMatrix.block(0,(_c+_channelOrder-1)*_nInputs,stackedChannelMatrix.rows(),(_e+1)*_nInputs),
                        stackedNoiseCovariance);
            } else
                softEstimations =  dynamic_cast <WithLinearDetectionParticleAddon *> (_particleFilter->getParticle(iParticle))->getLinearDetector(_estimatorIndex)->detect(stackedObservations,stackedChannelMatrix,stackedNoiseCovariance);

            // the evaluated proposal function (necessary for computing the weights) is initialized
            proposal = 1.0;

            // sampling
            for(iSampledSymbol=0;iSampledSymbol<(_nInputs*(_d+1));iSampledSymbol++)
            {
                s2q = dynamic_cast <WithLinearDetectionParticleAddon *> (_particleFilter->getParticle(iParticle))->getLinearDetector(_estimatorIndex)->nthSymbolVariance(iSampledSymbol,noiseVariances[iObservationToBeProcessed]);

                sumProb = 0.0;

                // the probability for each posible symbol alphabet is computed
                for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                {
                    symbolProb(iSampledSymbol,iAlphabet) = StatUtil::normalPdf(softEstimations(iSampledSymbol),dynamic_cast <WithLinearDetectionParticleAddon *> (_particleFilter->getParticle(iParticle))->getLinearDetector(_estimatorIndex)->nthSymbolGain(iSampledSymbol)*_alphabet[iAlphabet],s2q);

                    // the computed pdf is accumulated for normalizing purposes
                    sumProb += symbolProb(iSampledSymbol,iAlphabet);
                }

                if(sumProb!=0)
                    for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                    {
                        symbolProb(iSampledSymbol,iAlphabet) /= sumProb;
                    }
                else
                {
                    cout << "The sum of the probabilities is null." << endl;
                    for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                        symbolProb(iSampledSymbol,iAlphabet) = 1.0/double(_alphabet.length());
                }

                iSampled = StatUtil::discrete_rnd(symbolProb.row(iSampledSymbol));
                sampledSmoothingVector(iSampledSymbol) = _alphabet[iSampled];

                proposal *= symbolProb(iSampledSymbol,iSampled);
            }

            // sampled symbol vector is stored for the corresponding particle
            _particleFilter->getParticle(iParticle)->setSymbolVector(iObservationToBeProcessed,sampledSmoothingVector.start(_nInputs));

            // all the symbol vectors involved in the smoothing are kept in "forWeightUpdateNeededSymbols"
            // i) the already known:
            forWeightUpdateNeededSymbols.block(0,0,_nInputs,_channelOrder-1) = _particleFilter->getParticle(iParticle)->getSymbolVectors(iObservationToBeProcessed-_channelOrder+1,iObservationToBeProcessed-1);

            // ii) the just sampled
            forWeightUpdateNeededSymbols.block(0,_channelOrder-1,_nInputs,_d+1) = Util::toMatrix(sampledSmoothingVector,columnwise,_nInputs);

            likelihoodsProd = smoothedLikelihood(matricesToStack,forWeightUpdateNeededSymbols,iObservationToBeProcessed,observations,noiseVariances);

            // the weight is updated
            _particleFilter->getParticle(iParticle)->setWeight((likelihoodsProd/proposal)*_particleFilter->getParticle(iParticle)->getWeight());

            // and the estimation of the channel matrix
            dynamic_cast <WithChannelEstimationParticleAddon *> (_particleFilter->getParticle(iParticle))->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,
                                                dynamic_cast <WithChannelEstimationParticleAddon *> (_particleFilter->getParticle(iParticle))->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),
                                                    forWeightUpdateNeededSymbols.block(0,0,_nInputs,_channelOrder),noiseVariances[iObservationToBeProcessed]));            
        } // for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)

        _particleFilter->normalizeWeights();

        // if it's not the last time instant
        if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
            _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

void LinearFilterBasedSMCAlgorithm::beforeInitializingParticles(const MatrixXd &observations, const MatrixXd &trainingSequence)
{
    _linearDetector->stateStepsFromObservationsSequence(observations,_d,_preamble.cols(),_preamble.cols()+trainingSequence.cols());
}
