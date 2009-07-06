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

CDMAunknownActiveUsersSISopt::CDMAunknownActiveUsersSISopt(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances)
{
}


CDMAunknownActiveUsersSISopt::~CDMAunknownActiveUsersSISopt()
{
}


tMatrix CDMAunknownActiveUsersSISopt::getDetectedSymbolVectors()
{
    return SMCAlgorithm::getDetectedSymbolVectors();
}

vector< tMatrix > CDMAunknownActiveUsersSISopt::GetEstimatedChannelMatrices()
{
    return SMCAlgorithm::GetEstimatedChannelMatrices();
}

void CDMAunknownActiveUsersSISopt::BeforeInitializingParticles(const tMatrix& observations, const tMatrix& trainingSequence)
{
    SMCAlgorithm::BeforeInitializingParticles(observations, trainingSequence);
}

void CDMAunknownActiveUsersSISopt::InitializeParticles()
{
    SMCAlgorithm::InitializeParticles();
}

void CDMAunknownActiveUsersSISopt::Process(const tMatrix& observations, vector< double > noiseVariances)
{
    cout << "in Process: exiting..." << endl;
    exit(0);
    
    int k,iParticle,iSampledVector;
    vector<tSymbol> testedVector(_nInputs),sampledVector(_nInputs);
    tRange mMinus1FirstColumns(0,_channelOrder-2);

    // it selects all rows in the symbols Matrix
    tRange rAll;

    // it includes all symbol vectors involved in the smoothing
    tMatrix involvedSymbolVectors(_nInputs,_channelOrder);

    uint nSymbolVectors = (int) pow((double)_alphabet.length(),(double)_nInputs);

    // a likelihood is computed for every possible symbol vector
    vector<double> likelihoods(nSymbolVectors);

    tRange mPrecedentColumns(_startDetectionTime-_channelOrder+1,_startDetectionTime);
    tRange mMinus1PrecedentColumns(_startDetectionTime-_channelOrder+1,_startDetectionTime-1);

    double likelihoodsSum;

    // for each time instant
    for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
        for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)
        {
            ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

            // the m-1 already detected symbol vectors are copied into the matrix:
            involvedSymbolVectors(rAll,mMinus1FirstColumns).inject(processedParticle->GetSymbolVectors(mMinus1PrecedentColumns));

            for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
            {
                // the corresponding testing vector is generated from the index
                _alphabet.int2symbolsArray(iTestedVector,testedVector);

                // current tested vector is copied in the m-th position
                for(k=0;k<_nInputs;k++)
                    involvedSymbolVectors(k,_channelOrder-1) = testedVector[k];

                likelihoods[iTestedVector] = processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->likelihood(observations.col(iObservationToBeProcessed),involvedSymbolVectors,noiseVariances[iObservationToBeProcessed]);
            } // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

            likelihoodsSum = Util::sum(likelihoods);
            vector<double> probabilities = likelihoods;
            if(likelihoodsSum>0)
                // probabilities are computed by normalizing the likelihoods
                Util::normalize(probabilities);
            else
                // if all the likelihoods are null
                probabilities = vector<double>(nSymbolVectors,1.0/(double)nSymbolVectors);

            // one sample from the discrete distribution is taken
            iSampledVector = StatUtil::discrete_rnd(probabilities);

            // the above index is turned into a vector
            _alphabet.int2symbolsArray(iSampledVector,sampledVector);

            // sampled symbols are copied into the corresponding particle
            processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledVector);

            // channel matrix is estimated by means of the particle channel estimator
            processedParticle->SetChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),processedParticle->GetSymbolVectors(mPrecedentColumns),noiseVariances[iObservationToBeProcessed]));

            processedParticle->SetWeight(processedParticle->GetWeight()* Util::sum(likelihoods));
        } // for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)

        mPrecedentColumns = mPrecedentColumns + 1;
        mMinus1PrecedentColumns = mMinus1PrecedentColumns + 1;

        _particleFilter->NormalizeWeights();

        // if it's not the last time instant
        if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
//          Resampling();
            _resamplingAlgorithm->ResampleWhenNecessary(_particleFilter);

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

