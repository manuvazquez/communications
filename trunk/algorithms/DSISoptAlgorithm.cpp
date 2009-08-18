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
#include "DSISoptAlgorithm.h"

// #define DEBUG
// #define PRINT_INFO

DSISoptAlgorithm::DSISoptAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m, ChannelMatrixEstimator *channelEstimator, tMatrix preamble, int smoothingLag, int nParticles,ResamplingAlgorithm *resamplingAlgorithm, const tMatrix &channelMatrixMean, const tMatrix &channelMatrixVariances): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected,m,  channelEstimator, preamble, smoothingLag, nParticles,resamplingAlgorithm,channelMatrixMean,channelMatrixVariances)
{
//     _randomParticlesInitilization = true;
}

// void DSISoptAlgorithm::process(const tMatrix &observations, vector< double > noiseVariances)
void DSISoptAlgorithm::process(const MatrixXd &observations, vector< double > noiseVariances)
{
	uint k,iSmoothingVector;
	int iSmoothingLag,iParticle,iSampledVector;
	vector<tSymbol> testedVector(_nInputs),testedSmoothingVector(_nInputs*_d),sampledVector(_nInputs);
	double auxLikelihoodsProd;
	ChannelMatrixEstimator *channelEstimatorClone;

	// it selects all rows in the symbols Matrix
	tRange rAll;

	// it includes all symbol vectors involved in the smoothing
    MatrixXd smoothingSymbolVectors(_nInputs,_channelOrder+_d);

	uint nSymbolVectors = (int) pow((double)_alphabet.length(),(double)_nInputs);
	uint nSmoothingVectors = (int) pow((double)_alphabet.length(),(double)(_nInputs*_d));

	// a likelihood is computed for every possible symbol vector
    VectorXd likelihoods(nSymbolVectors);

	// for each time instant
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
	{
		for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
		{
#ifdef PRINT_INFO
            cout << "iObservationToBeProcessed = " << iObservationToBeProcessed << " iParticle = " << iParticle << endl;
#endif      
			ParticleWithChannelEstimation *processedParticle = dynamic_cast<ParticleWithChannelEstimation *>(_particleFilter->getParticle(iParticle));

			// the m-1 already detected symbol vectors are copied into the matrix:
            smoothingSymbolVectors.block(0,0,_nInputs,_channelOrder-1) = processedParticle->getSymbolVectors().block(0,iObservationToBeProcessed-_channelOrder+1,_nInputs,_channelOrder-1);

			for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.int2symbolsArray(iTestedVector,testedVector);

				// current tested vector is copied in the m-th position
				for(k=0;k<static_cast<uint>(_nInputs);k++)
					smoothingSymbolVectors(k,_channelOrder-1) = testedVector[k];

				likelihoods(iTestedVector) = 0.0;

				// every possible smoothing sequence is tested
				for(iSmoothingVector=0;iSmoothingVector<nSmoothingVectors;iSmoothingVector++)
				{
					// a testing smoothing vector is generated for the index
					_alphabet.int2symbolsArray(iSmoothingVector,testedSmoothingVector);

					// symbols used for smoothing are copied into "smoothingSymbolVectors"
					for(k=0;k<testedSmoothingVector.size();k++)
						smoothingSymbolVectors((_nInputsXchannelOrder+k)%_nInputs,(_nInputsXchannelOrder+k)/_nInputs) = testedSmoothingVector[k];

					// a clone of the channel estimator is generated because this must not be modified
					channelEstimatorClone = processedParticle->getChannelMatrixEstimator(_estimatorIndex)->clone();

					auxLikelihoodsProd = 1.0;

					for(iSmoothingLag=0;iSmoothingLag<=_d;iSmoothingLag++)
					{

						// the likelihood is computed and accumulated
                        auxLikelihoodsProd *= channelEstimatorClone->likelihood(observations.col(iObservationToBeProcessed+iSmoothingLag),smoothingSymbolVectors.block(0,iSmoothingLag,_nInputs,_channelOrder),noiseVariances[iObservationToBeProcessed+iSmoothingLag]);

						// a step in the Kalman Filter
                        channelEstimatorClone->nextMatrix(observations.col(iObservationToBeProcessed+iSmoothingLag),smoothingSymbolVectors.block(0,iSmoothingLag,_nInputs,_channelOrder),noiseVariances[iObservationToBeProcessed+iSmoothingLag]);

					} // for(iSmoothingLag=0;iSmoothingLag<=_d;iSmoothingLag++)

					// memory of the clone is freed
					delete channelEstimatorClone;

					// the likelihood is accumulated in the proper position of vector:
					likelihoods(iTestedVector) += auxLikelihoodsProd;
				} // for(iSmoothingVector=0;iSmoothingVector<nSmoothingVectors;iSmoothingVector++)

			} // for(int iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

			VectorXd probabilities(nSymbolVectors);
			try {
				// probabilities are computed by normalizing the likelihoods
				probabilities = Util::normalize(likelihoods);
			} catch (AllElementsNullException) {
				// if all the likelihoods are null
                probabilities.setConstant(1.0/(double)nSymbolVectors);
			}

			// one sample from the discrete distribution is taken
			iSampledVector = StatUtil::discrete_rnd(probabilities);

			// the above index is turned into a vector
			_alphabet.int2symbolsArray(iSampledVector,sampledVector);

			// sampled symbols are copied into the corresponding particle
			processedParticle->setSymbolVector(iObservationToBeProcessed,sampledVector);

			// channel matrix is estimated by means of the particle channel estimator
            VectorXd currentObservation = observations.col(iObservationToBeProcessed);
            MatrixXd currentSymbolVectors = processedParticle->getSymbolVectors().block(0,iObservationToBeProcessed-_channelOrder+1,_nInputs,_channelOrder);
            processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(Util::eigen2lapack(currentObservation),Util::eigen2lapack(currentSymbolVectors),noiseVariances[iObservationToBeProcessed]));

			processedParticle->setWeight(processedParticle->getWeight()* Util::sum(likelihoods));

		} // for(iParticle=0;iParticle<_nParticles;iParticle++)

		_particleFilter->normalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
            _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);
	} // for each time instant
}

