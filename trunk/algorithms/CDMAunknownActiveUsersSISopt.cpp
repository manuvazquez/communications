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
#include "CDMAunknownActiveUsersSISopt.h"

// #define IMPORT_REAL_DATA

#ifdef IMPORT_REAL_DATA
    extern MIMOChannel *realChannel;
    extern tMatrix *realSymbols;
    extern Noise *realNoise;
#endif

// #define DEBUG_MSE_THRESHOLD 0.5
// #define DEBUG
// #define DEBUG_CHANNEL_SAMPLES

// #define DEBUG5

CDMAunknownActiveUsersSISopt::CDMAunknownActiveUsersSISopt(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances,const double userPersistenceProb,const double newActiveUserProb,const double userPriorProb): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances),_userPersistenceProb(userPersistenceProb),_newActiveUserProb(newActiveUserProb),_userPriorProb(userPriorProb)
{    
    _randomParticlesInitilization = true;    
}

// vector< tMatrix > CDMAunknownActiveUsersSISopt::getEstimatedChannelMatrices()
// {
//     return SMCAlgorithm::getEstimatedChannelMatrices();
// }

void CDMAunknownActiveUsersSISopt::InitializeParticles()
{
    ChannelMatrixEstimator *channelMatrixEstimatorClone;
    tVector channelMean = Util::toVector(_channelMatrixMean,rowwise);
    tMatrix channelCovariance = LaGenMatDouble::from_diag(_channelMatrixVariances);

    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
    {
        channelMatrixEstimatorClone = _channelEstimator->clone();
        
        if(_randomParticlesInitilization)
            channelMatrixEstimatorClone->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::RandMatrix(channelMean,channelCovariance),rowwise,_Nr));
        
        _particleFilter->addParticle(new ParticleWithChannelEstimationAndActiveUsers(1.0/(double)_particleFilter->capacity(),_nInputs,_iLastSymbolVectorToBeDetected,channelMatrixEstimatorClone));

        // if there is preamble...
        if(_preamble.cols()!=0)
            _particleFilter->getParticle(iParticle)->setSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);
    }
}

void CDMAunknownActiveUsersSISopt::Process(const tMatrix& observations, vector< double > noiseVariances)
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

    // it selects all rows in the symbols Matrix
    tRange rAll;

    // tVector containing the symbols
    tVector symbolsVector(_nInputs);

    // a likelihood is computed for every possible symbol vector
    vector<double> testedCombination(_nInputs),likelihoods(nCombinations);
    
    // for storing users activity
    vector<bool> usersActivity(_nInputs);

    double likelihoodsSum;

    // for each time instant
    for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
#ifdef DEBUG
        cout << "******************************************** iObservationToBeProcessed = " << iObservationToBeProcessed << " **************************" << endl;
#endif    
        for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
        {
#ifdef DEBUG
            cout << "------------- iObservationToBeProcessed = " << iObservationToBeProcessed << " (iParticle = " << iParticle << ") -------------" << endl;
#endif        
            ParticleWithChannelEstimationAndActiveUsers *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndActiveUsers *> (_particleFilter->getParticle(iParticle));

#ifdef DEBUG_CHANNEL_SAMPLES
            tMatrix channelMatrixEstimation = processedParticle->getChannelMatrixEstimator(_estimatorIndex)->lastEstimatedChannelMatrix();
            cout << "channel matrix " << (*realChannel)[iObservationToBeProcessed];
            cout << "last estimated channel matrix " << channelMatrixEstimation;
            double normalizedMSE = Util::normalizedSquareError(channelMatrixEstimation,(*realChannel)[iObservationToBeProcessed]);
            cout << "normalized MSE = " << normalizedMSE << endl;
            if(normalizedMSE < DEBUG_MSE_THRESHOLD)
                getchar();            
#endif

            for(uint iTestedCombination=0;iTestedCombination<nCombinations;iTestedCombination++)
            {
                // the corresponding testing vector is generated from the index
                extendedAlphabet.int2symbolsArray(iTestedCombination,testedCombination);

                // symbols are copied into a lapackpp vector
                for(k=0;k<_nInputs;k++)
                    symbolsVector(k) = testedCombination[k];

                likelihoods[iTestedCombination] = processedParticle->getChannelMatrixEstimator(_estimatorIndex)->likelihood(observations.col(iObservationToBeProcessed),symbolsVector,noiseVariances[iObservationToBeProcessed]);

#ifdef DEBUG5
                cout << "probability of users active times symbol vector = " << probSymbolsVectorXprobActiveUsers(symbolsVector) << endl;
#endif

                // the probability of these users being active and this particular symbol vector being transmitted is computed...
                // ...either taking into account the users previous state in case this exists
                if(iObservationToBeProcessed!=_startDetectionTime)
                    likelihoods[iTestedCombination] *=probSymbolsVectorXprobActiveUsers(symbolsVector,processedParticle->getActivityAtTime(iObservationToBeProcessed-1));
                // ...or not doing so when it's the first time instant
                else
                    likelihoods[iTestedCombination] *=probSymbolsVectorXprobActiveUsers(symbolsVector);
                
            } // for(uint iTestedCombination=0;iTestedCombination<nCombinations;iTestedCombination++)

            likelihoodsSum = Util::sum(likelihoods);
            vector<double> probabilities = likelihoods;
            if(likelihoodsSum>0)
                // probabilities are computed by normalizing the likelihoods
                Util::normalize(probabilities);
            else
                // if all the likelihoods are null
                probabilities = vector<double>(nCombinations,1.0/(double)nCombinations);

#ifdef DEBUG
            if(iParticle>=0)
            {
                cout << "computed probabilities" << endl;
                cout << "observations" << endl << observations.col(iObservationToBeProcessed);
                cout << "noise affecting..." << endl << (*realNoise)[iObservationToBeProcessed];
                cout << "true transmitted symbol vector" << endl << realSymbols->col(iObservationToBeProcessed);
                for(uint iTestedCombination=0;iTestedCombination<nCombinations;iTestedCombination++)
                {
                    extendedAlphabet.int2symbolsArray(iTestedCombination,testedCombination);
                    cout << "tested combination is   ";
                    Util::print(testedCombination);
                    cout << " its probability = " << probabilities[iTestedCombination] << endl;
                }
                if(normalizedMSE < DEBUG_MSE_THRESHOLD)
                    getchar();         
            }
#endif

            // one sample from the discrete distribution is taken
            iSampledVector = StatUtil::discrete_rnd(probabilities);
            
//             iSampledVector = Util::max(probabilities); // <---------------------------------------------------------------------------------

            // the above index is turned into a vector
            extendedAlphabet.int2symbolsArray(iSampledVector,sampledVector);
                   
            // sampled symbols are copied into the corresponding particle
            processedParticle->setSymbolVector(iObservationToBeProcessed,sampledVector);

            // channel matrix is estimated by means of the particle channel estimator
            processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),processedParticle->getSymbolVector(iObservationToBeProcessed),noiseVariances[iObservationToBeProcessed]));
                        
#ifdef DEBUG
            cout << "sampled symbol vector by particle " << iParticle << " ";
            Util::print(sampledVector);            
            cout << endl;
            cout << "updated channel matrix" << endl << processedParticle->getChannelMatrixEstimator(_estimatorIndex)->lastEstimatedChannelMatrix();
            if(normalizedMSE < DEBUG_MSE_THRESHOLD)
                getchar();             
#endif                             
                        
            // users activity
            for(k=0;k<_nInputs;k++)
                usersActivity[k] = isUserActive(sampledVector[k]);
            
            processedParticle->setActivityAtTime(iObservationToBeProcessed,usersActivity);
            
            processedParticle->setWeight(processedParticle->getWeight()* Util::sum(likelihoods));
        } // for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)

        _particleFilter->normalizeWeights();

#ifdef DEBUG
//         cout << "weights" << endl << _particleFilter->getWeightsVector();
        _particleFilter->printWeights();
        getchar();
#endif

        // if it's not the last time instant
        if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
            _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);        
//             cout << "resampled = " << _resamplingAlgorithm->resampleWhenNecessary(_particleFilter) << endl;

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

// double CDMAunknownActiveUsersSISopt::probSymbolsVectorGivenActiveUsers(const tVector &v) const
// {
//     uint nActiveUsers = 0;
//     
//     for(uint i=0;i<v.size();i++)
//         nActiveUsers = nActiveUsers + (v(i)==0.0);
// 
//     return 1.0/pow((double)_alphabet.length(),(double)nActiveUsers);
// }

double CDMAunknownActiveUsersSISopt::probSymbolsVectorXprobActiveUsers(const tVector &symbolsVector, const std::vector<bool> &lastUsersActivity) const
{
    if(symbolsVector.size()!=_nInputs)
        throw RuntimeException("CDMAunknownActiveUsersSISopt::probSymbolsVectorXprobActiveUsers: symbols vector dimensions are wrong.");
    
    if(static_cast<uint> (symbolsVector.size())!=lastUsersActivity.size())
        throw RuntimeException("CDMAunknownActiveUsersSISopt::probSymbolsVectorXprobActiveUsers: symbols vector size doesn't coincide with that of the vector containing information about the users activity in the previous time instant.");
        
    double probUsersActivity = 1.0;
    double probSymbolsVector = 1.0;
    
    for(int i=0;i<symbolsVector.size();i++)
    {
        if(isUserActive(symbolsVector(i)))
        {
            probUsersActivity *= lastUsersActivity[i]?_userPersistenceProb:_newActiveUserProb;
            probSymbolsVector /= double(_alphabet.length());
        }
        else
            probUsersActivity *= lastUsersActivity[i]?(1.0-_userPersistenceProb):(1.0-_newActiveUserProb);
    }

    return probSymbolsVector*probUsersActivity;
}

double CDMAunknownActiveUsersSISopt::probSymbolsVectorXprobActiveUsers(const tVector &symbolsVector) const
{
    if(symbolsVector.size()!=_nInputs)
        throw RuntimeException("CDMAunknownActiveUsersSISopt::probSymbolsVectorXprobActiveUsers: symbols vector dimensions are wrong.");
                        
    double probUsersActivity = 1.0;
    double probSymbolsVector = 1.0;
    
    for(int i=0;i<symbolsVector.size();i++)
    {
        if(isUserActive(symbolsVector(i)))
        {
            probUsersActivity *= _userPriorProb;
            probSymbolsVector /= double(_alphabet.length());
        }
        else
            probUsersActivity *= 1 - _userPriorProb;
    }

    return probSymbolsVector*probUsersActivity;
}
