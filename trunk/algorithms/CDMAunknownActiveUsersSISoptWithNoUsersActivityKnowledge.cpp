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
#include "CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge.h"

// #define IMPORT_REAL_DATA

#ifdef IMPORT_REAL_DATA
    extern MIMOChannel *realChannel;
    extern MatrixXd *realSymbols;
    extern Noise *realNoise;
#endif

CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge::CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const MatrixXd& channelMatrixMean, const MatrixXd& channelMatrixVariances): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances)
{    
    _randomParticlesInitilization = true;    
}


void CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge::initializeParticles()
{
    ChannelMatrixEstimator *channelMatrixEstimatorClone;

    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
    {
        channelMatrixEstimatorClone = _channelEstimator->clone();
        
        if(_randomParticlesInitilization)
            channelMatrixEstimatorClone->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::randnMatrix(_channelMean,_channelCovariance),rowwise,_Nr));
        
        _particleFilter->addParticle(new ParticleWithChannelEstimation(1.0/(double)_particleFilter->capacity(),_nInputs,_iLastSymbolVectorToBeDetected,channelMatrixEstimatorClone));

        // if there is preamble...
        if(_preamble.cols()!=0)
            _particleFilter->getParticle(iParticle)->setSymbolVectors(0,_preamble.cols(),_preamble);
    }
}

void CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge::process(const MatrixXd& observations, vector< double > noiseVariances)
{    
    // a new alphabet extended with 0 (that meaning, no symbol is transmitted)
    vector<tSymbol> extendedAlphabetSymbols(_alphabet.length()+1);
    
    for(int i=0;i<_alphabet.length();i++)
        extendedAlphabetSymbols[i] = _alphabet[i];
    extendedAlphabetSymbols[_alphabet.length()] = 0.0;
    
    Alphabet extendedAlphabet(extendedAlphabetSymbols);
        
//     extendedAlphabet = _alphabet; // <-----------------------------------------------------------
    
    uint nCombinations = (int) pow((double)(extendedAlphabet.length()),(double)_nInputs);
    
    vector<tSymbol> combination(_nInputs,_alphabet[0]);

    int k,iParticle,iSampledVector;
    vector<tSymbol> sampledVector(_nInputs);

    // tVector containing the symbols
    VectorXd symbolsVector(_nInputs);

    // a likelihood is computed for every possible symbol vector
    vector<double> testedCombination(_nInputs),likelihoods(nCombinations);

    double likelihoodsSum;

    // for each time instant
    for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
        for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
        {
            ParticleWithChannelEstimation *processedParticle = dynamic_cast<ParticleWithChannelEstimation *> (_particleFilter->getParticle(iParticle));

            for(uint iTestedCombination=0;iTestedCombination<nCombinations;iTestedCombination++)
            {
                // the corresponding testing vector is generated from the index
                extendedAlphabet.int2symbolsArray(iTestedCombination,testedCombination);

                // symbols are copied into a vector
                for(k=0;k<_nInputs;k++)
                    symbolsVector(k) = testedCombination[k];

                likelihoods[iTestedCombination] = processedParticle->getChannelMatrixEstimator(_estimatorIndex)->likelihood(observations.col(iObservationToBeProcessed),symbolsVector,noiseVariances[iObservationToBeProcessed]);
                
            } // for(uint iTestedCombination=0;iTestedCombination<nCombinations;iTestedCombination++)

            likelihoodsSum = Util::sum(likelihoods);
            vector<double> probabilities = likelihoods;
            if(likelihoodsSum>0)
                // probabilities are computed by normalizing the likelihoods
                Util::normalize(probabilities);
            else
                // if all the likelihoods are null
                probabilities = vector<double>(nCombinations,1.0/(double)nCombinations);

            // one sample from the discrete distribution is taken
            iSampledVector = StatUtil::discrete_rnd(probabilities);
            
//             iSampledVector = Util::max(probabilities); // <---------------------------------------------------------------------------------

            // the above index is turned into a vector
            extendedAlphabet.int2symbolsArray(iSampledVector,sampledVector);
                   
            // sampled symbols are copied into the corresponding particle
            processedParticle->setSymbolVector(iObservationToBeProcessed,sampledVector);

            // channel matrix is estimated by means of the particle channel estimator
            processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),processedParticle->getSymbolVector(iObservationToBeProcessed),noiseVariances[iObservationToBeProcessed]));
            
            processedParticle->setWeight(processedParticle->getWeight()* likelihoodsSum);
        } // for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)

        _particleFilter->normalizeWeights();

        // if it's not the last time instant
        if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
            _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);        

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}
