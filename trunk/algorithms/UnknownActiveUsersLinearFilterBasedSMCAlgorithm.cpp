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

// #define DEBUG
// #define DEBUG_CHANNEL_SAMPLES
// #define DEBUG_MSE_THRESHOLD 1.5

#define IMPORT_REAL_DATA

#ifdef IMPORT_REAL_DATA
    extern MIMOChannel *realChannel;
    extern tMatrix *realSymbols;
    extern Noise *realNoise;
#endif

UnknownActiveUsersLinearFilterBasedSMCAlgorithm::UnknownActiveUsersLinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet, int L, int Nr, int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, LinearDetector *linearDetector, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances, const UsersActivityDistribution &usersActivityPdf): SMCAlgorithm(name, alphabet, L, Nr, N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances),_linearDetector(linearDetector->clone()),_usersActivityPdf(usersActivityPdf)
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
            channelMatrixEstimatorClone->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::randMatrix(channelMean,channelCovariance),rowwise,_Nr));
        
        _particleFilter->addParticle(new ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers(1.0/(double)_particleFilter->capacity(),_nInputs,_iLastSymbolVectorToBeDetected,channelMatrixEstimatorClone,_linearDetector->clone()));

        if(_preamble.cols()!=0)
            _particleFilter->getParticle(iParticle)->setSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);
    }
}

void UnknownActiveUsersLinearFilterBasedSMCAlgorithm::process(const tMatrix& observations, vector< double > noiseVariances)
{
    int iParticle,iSampledSymbol,iAlphabet,iSampled;
    tRange rAll;
    double proposal,s2q,sumProb,likelihood;
    tMatrix forWeightUpdateNeededSymbols(_nInputs,_channelOrder+_d);
    
    int iOutput,iInput;
    tMatrix channelMatrixSample;
    tVector sampledVector(_nInputs);
    tMatrix symbolProb(_nInputs,_alphabet.length());
    tVector predictedNoiselessObservation(_nOutputs);
    double probSymbolsVectorApriori;

    for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
        // noise covariance needs to be constructed
        tMatrix noiseCovariance = LaGenMatDouble::zeros(_nOutputs,_nOutputs);
        for(iOutput=0;iOutput<_nOutputs;iOutput++)
            noiseCovariance(iOutput,iOutput) = noiseVariances[iObservationToBeProcessed];

        for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
        {
#ifdef DEBUG
            cout << "------ iObservationToBeProcessed = " << iObservationToBeProcessed << " iParticle = " << iParticle << " -----------" << endl;
#endif
            ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers *> (_particleFilter->getParticle(iParticle));            
            
            channelMatrixSample = (dynamic_cast<KalmanEstimator *> (_particleFilter->getParticle(iParticle)->getChannelMatrixEstimator(_estimatorIndex)))->sampleFromPredictive();            

#ifdef DEBUG_CHANNEL_SAMPLES
            tMatrix channelMatrixEstimation = processedParticle->getChannelMatrixEstimator(_estimatorIndex)->lastEstimatedChannelMatrix();
            cout << "channel matrix " << (*realChannel)[iObservationToBeProcessed];
            cout << "last estimated channel matrix " << channelMatrixEstimation;
            cout << "channelMatrixSample = " << endl << channelMatrixSample;
            double normalizedMSE = Util::normalizedSquareError(channelMatrixEstimation,(*realChannel)[iObservationToBeProcessed]);
            cout << "normalized MSE = " << normalizedMSE << endl;
            if(normalizedMSE < DEBUG_MSE_THRESHOLD)
                getchar();                    
#endif

            // sampling of the users activity using:            
            // i) a priori (first time instant)
            if(iObservationToBeProcessed==_startDetectionTime)
                for(iInput=0;iInput<_nInputs;iInput++)
                    processedParticle->setUserActivity(iInput,iObservationToBeProcessed,_usersActivityPdf.sampleFromPrior());
            // ii) conditional pdf (after first time instant)
            else
                for(iInput=0;iInput<_nInputs;iInput++)
                    processedParticle->setUserActivity(iInput,iObservationToBeProcessed,_usersActivityPdf.sampleGivenItWas(processedParticle->getUserActivity(iInput,iObservationToBeProcessed-1)));                

#ifdef DEBUG
            if(iObservationToBeProcessed>_startDetectionTime)
            {
                cout << "previous users activity" << endl;
                Util::print(processedParticle->getActivityAtTime(iObservationToBeProcessed-1));
                cout << endl;
            }
            cout << "sampled users activity" << endl;
            Util::print(processedParticle->getActivityAtTime(iObservationToBeProcessed));
            cout << endl;
#endif

            tVector softEstimations;

            // the sampled channel matrix is used to obtain soft estimations of the transmitted symbols
            softEstimations =  processedParticle->getLinearDetector(_estimatorIndex)->detect(observations.col(iObservationToBeProcessed),channelMatrixSample,noiseCovariance);

            // the evaluated proposal function (necessary for computing the weights) is initialized
            proposal = 1.0;
            
            // and so it is the a priori prob
            probSymbolsVectorApriori = 1.0;

            // sampling
            for(iSampledSymbol=0;iSampledSymbol<_nInputs;iSampledSymbol++)
            {
                // if the user is not active
                if(!processedParticle->getUserActivity(iSampledSymbol,iObservationToBeProcessed))
                {
                    // the symbol is known to be 0.0 (and the proposal doesn' have to be updated)
                    sampledVector(iSampledSymbol) = 0.0;
                    continue;
                }
                
                s2q = processedParticle->getLinearDetector(_estimatorIndex)->nthSymbolVariance(iSampledSymbol);

#ifdef DEBUG
                cout << "s2q = " << s2q << endl;
#endif

                sumProb = 0.0;

                // the probability for each posible symbol alphabet is computed
                for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                {
                    symbolProb(iSampledSymbol,iAlphabet) = StatUtil::normalPdf(softEstimations(iSampledSymbol),processedParticle->getLinearDetector(_estimatorIndex)->nthSymbolGain(iSampledSymbol)*_alphabet[iAlphabet],s2q);

                    // the computed pdf is accumulated for normalizing purposes
                    sumProb += symbolProb(iSampledSymbol,iAlphabet);
                }

#ifdef DEBUG
//                 if(iObservationToBeProcessed==5 && iParticle==3)
//                     cout << "sumProb = " << sumProb << endl;
//                 if (sumProb==0)
//                     cout << "is zero!!!" << endl;
#endif

//                 try {
                if(sumProb!=0)
                    for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                        symbolProb(iSampledSymbol,iAlphabet) /= sumProb;
                else
                {
//                 }catch(exception e){
                    cout << "The sum of the probabilities is null." << endl;
                    for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                        symbolProb(iSampledSymbol,iAlphabet) = 1.0/double(_alphabet.length());
                }

                iSampled = StatUtil::discrete_rnd(symbolProb.row(iSampledSymbol));
                sampledVector(iSampledSymbol) = _alphabet[iSampled];

                proposal *= symbolProb(iSampledSymbol,iSampled);
                probSymbolsVectorApriori /= double(_alphabet.length());
            }

#ifdef DEBUG
            cout << "soft estimations" << endl << softEstimations;
            cout << "true transmitted symbol vector" << endl << realSymbols->col(iObservationToBeProcessed);
            cout << "sampled vector" << endl << sampledVector;
            cout << "from" << endl << symbolProb;
            if(normalizedMSE < DEBUG_MSE_THRESHOLD)
                getchar();            
#endif

            // sampled symbol vector is stored for the corresponding particle
            processedParticle->setSymbolVector(iObservationToBeProcessed,sampledVector);

            // predictedNoiselessObservation = matricesToStack[iSmoothing] * stackedSymbolVector
            Blas_Mat_Vec_Mult(channelMatrixSample,sampledVector,predictedNoiselessObservation);

            likelihood = StatUtil::normalPdf(observations.col(iObservationToBeProcessed),predictedNoiselessObservation,noiseVariances[iObservationToBeProcessed]);

            // the weight is updated
            processedParticle->setWeight((likelihood*probSymbolsVectorApriori/proposal)*processedParticle->getWeight());

            // and the estimation of the channel matrix
            processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,
                                                processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),sampledVector,noiseVariances[iObservationToBeProcessed]));

#ifdef DEBUG
//             getchar();
#endif

        } // for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)

        _particleFilter->normalizeWeights();

        // if it's not the last time instant
        if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
            _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

