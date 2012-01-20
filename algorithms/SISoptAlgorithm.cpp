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
#include "SISoptAlgorithm.h"

SISoptAlgorithm::SISoptAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm, const MatrixXd& channelMatrixMean, const MatrixXd& channelMatrixVariances): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, 0, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances)
{
}

// eigen
void SISoptAlgorithm::process(const MatrixXd& observations, vector<double> noiseVariances)
{
    uint k,iParticle,iSampledVector;
    vector<tSymbol> testedVector(_nInputs),sampledVector(_nInputs);

    // it includes all symbol vectors involved in the smoothing
    MatrixXd involvedSymbolVectors(_nInputs,_channelOrder);

    uint nSymbolVectors = uint(pow((double)_alphabet.length(),(double)_nInputs));

    // a likelihood is computed for every possible symbol vector
    VectorXd likelihoods(nSymbolVectors);

    // for each time instant
    for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
        for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
        {
            ParticleWithChannelEstimation *processedParticle = dynamic_cast<ParticleWithChannelEstimation *>(_particleFilter->getParticle(iParticle));

            // the m-1 already detected symbol vectors are copied into the matrix:
            involvedSymbolVectors.block(0,0,_nInputs,_channelOrder-1) = processedParticle->getSymbolVectors(iObservationToBeProcessed-_channelOrder+1,iObservationToBeProcessed-1);

            for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
            {
                // the corresponding testing vector is generated from the index
                _alphabet.int2symbolsArray(iTestedVector,testedVector);

                // current tested vector is copied in the m-th position
                for(k=0;k<_nInputs;k++)
                    involvedSymbolVectors(k,_channelOrder-1) = testedVector[k];

                likelihoods(iTestedVector) = processedParticle->getChannelMatrixEstimator(_estimatorIndex)->likelihood(observations.col(iObservationToBeProcessed),involvedSymbolVectors,noiseVariances[iObservationToBeProcessed]);
            } // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

            VectorXd probabilities(nSymbolVectors);
            try {
                // probabilities are computed by normalizing the likelihoods
                probabilities = Util::normalize(likelihoods);
            } catch (AllElementsNullException) {
                // if all the likelihoods are null
                probabilities /= (double)nSymbolVectors;
            }

            // one sample from the discrete distribution is taken
            iSampledVector = StatUtil::discrete_rnd(probabilities);

            // the above index is turned into a vector
            _alphabet.int2symbolsArray(iSampledVector,sampledVector);

            // sampled symbols are copied into the corresponding particle
            processedParticle->setSymbolVector(iObservationToBeProcessed,sampledVector);

            // channel matrix is estimated by means of the particle channel estimator
            processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),processedParticle->getSymbolVectors(iObservationToBeProcessed-_channelOrder+1,iObservationToBeProcessed),noiseVariances[iObservationToBeProcessed]));

//             processedParticle->setWeight(processedParticle->getWeight()* Util::sum(likelihoods));
			processedParticle->setWeight(processedParticle->getWeight()*likelihoods.sum());
        } // for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)

        _particleFilter->normalizeWeights();

        // if it's not the last time instant
        if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
            _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);

    } // for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}
