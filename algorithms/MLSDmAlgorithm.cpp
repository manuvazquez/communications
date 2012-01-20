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
#include "MLSDmAlgorithm.h"

// #define DEBUG

MLSDmAlgorithm::MLSDmAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, MatrixXd preamble, uint iFirstObservation, uint smoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance): MultipleChannelEstimatorsPerParticleSMCAlgorithm (name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_particleFilter(new ParticleFilter(nParticles)),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_particlesBestChannelOrders(nParticles)
{
}

MLSDmAlgorithm::~MLSDmAlgorithm()
{
    delete _particleFilter;
}

void MLSDmAlgorithm::initializeParticles()
{
    vector<ChannelMatrixEstimator *> channelEstimatorsClone(_channelEstimators.size());
    for(uint i=0;i<_candidateOrders.size();i++)
        channelEstimatorsClone[i] = _channelEstimators[i]->clone();

    // we begin with only one particle
    ParticleWithChannelEstimationAndChannelOrderAPP *particle = new ParticleWithChannelEstimationAndChannelOrderAPP(1.0,_nInputs,_iLastSymbolVectorToBeDetected+_d,channelEstimatorsClone);

//     particle->setSymbolVectors(0,_preamble.cols(),Util::lapack2eigen(_preamble));
    particle->setSymbolVectors(0,_preamble.cols(),_preamble);

    // the available APP's just before the _startDetectionTime instant are copied into the particle
    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        particle->setChannelOrderAPP(_channelOrderAPPs[0](iChannelOrder,_startDetectionTime-1),iChannelOrder);

    _particleFilter->addParticle(particle);
}

void MLSDmAlgorithm::process(const MatrixXd& observations, vector<double> noiseVariances)
{
    uint nSymbolVectors = (uint) pow((double)_alphabet.length(),(double)_nInputs);
    vector<tSymbol> testedVector(_nInputs);
    VectorXd computedObservations(_nOutputs);
    uint iCandidate,m,iBestUnnormalizedChannelOrderAPP,k,iParticle;
    uint iChannelOrder,iTestedVector;
    ParticleWithChannelEstimationAndChannelOrderAPP *processedParticle;
	
	// over all the survivor candidates
    double normConst;
	
	double likelihood;

    typedef struct{
        uint fromParticle;
        MatrixXd symbolVectorsMatrix;
        uint iBestChannelOrder;
        VectorXd unnormalizedChannelOrderAPPs;
        double likelihood;
        double weight;
    }tParticleCandidate;

    tParticleCandidate *particleCandidates = new tParticleCandidate[_particleFilter->capacity()*nSymbolVectors];

    // "symbolVectorsMatrix" will contain all the symbols involved in the current observation
    MatrixXd symbolVectorsMatrix(_nInputs,_maxOrder);
    VectorXd symbolsVector;

    uint lastSymbolVectorStart = _nInputsXmaxChannelOrder - _nInputs;

    vector<bool> activeCandidateOrders(_candidateOrders.size(),true);
    uint iBestChannelOrder = 0,timesBestChannelOrder = 0;

    for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected+_d;iObservationToBeProcessed++)
    {
        // it keeps track of the place where a new tParticleCandidate will be stored within the array
        iCandidate = 0;

        normConst = 0.0;

        // the candidates from all the particles are generated
        for(iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)
        {
            processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->getParticle(iParticle));

            symbolVectorsMatrix.block(0,0,_nInputs,_maxOrder-1) = processedParticle->getSymbolVectors(iObservationToBeProcessed-_maxOrder+1,iObservationToBeProcessed-1);
            symbolsVector = Util::toVector(symbolVectorsMatrix,columnwise);

            for(iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
            {
                // the corresponding testing vector is generated from the index
                _alphabet.int2symbolsArray(iTestedVector,testedVector);

                // current tested vector is copied in the m-th position
                for(k=0;k<_nInputs;k++)
                    symbolVectorsMatrix(k,_maxOrder-1) = symbolsVector(lastSymbolVectorStart+k) = testedVector[k];

                likelihood = 0.0;

                particleCandidates[iCandidate].unnormalizedChannelOrderAPPs = VectorXd(_candidateOrders.size());

                for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
                {
                    m = _candidateOrders[iChannelOrder];

                    MatrixXd involvedSymbolVectors = symbolVectorsMatrix.block(0,_maxOrder-m,_nInputs,m);

                    // the AR coefficiente is accounted for
//                  involvedSymbolVectors *= _ARcoefficient; <--------------------------------- (when Kalman estimator, it is already accounted for)

                    particleCandidates[iCandidate].unnormalizedChannelOrderAPPs(iChannelOrder) = processedParticle->getChannelOrderAPP(iChannelOrder)*
								  processedParticle->getChannelMatrixEstimator(iChannelOrder)->likelihood(observations.col(iObservationToBeProcessed),involvedSymbolVectors,noiseVariances[iObservationToBeProcessed]);


                    likelihood += particleCandidates[iCandidate].unnormalizedChannelOrderAPPs(iChannelOrder);
                }

                // if the likelihood is zero, we don't generate a candidate for this particle and this symbol vector
                if(likelihood==0.0)
                    continue;

                particleCandidates[iCandidate].unnormalizedChannelOrderAPPs.maxCoeff(&iBestUnnormalizedChannelOrderAPP);
                particleCandidates[iCandidate].fromParticle = iParticle;
                particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;
                particleCandidates[iCandidate].iBestChannelOrder = iBestUnnormalizedChannelOrderAPP;
                particleCandidates[iCandidate].likelihood = likelihood;
                particleCandidates[iCandidate].weight = processedParticle->getWeight()*likelihood;
                normConst += particleCandidates[iCandidate].weight;

                iCandidate++;
            } // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

        } // for(int iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)

        // if none of the candidates was valid
        if(iCandidate==0)
        {
            VectorXd uniformDistribution(_alphabet.length());
            for(uint iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                uniformDistribution(iAlphabet) = 1.0/_alphabet.length();

            for(iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)
            {
                processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->getParticle(iParticle));
                
                symbolVectorsMatrix.block(0,0,_nInputs,_maxOrder-1) = processedParticle->getSymbolVectors(iObservationToBeProcessed-_maxOrder+1,iObservationToBeProcessed-1);
                
				// the transmitted symbols are sampled from a uniform distribution
                for(k=0;k<_nInputs;k++)
                    symbolVectorsMatrix(k,_maxOrder-1) = _alphabet[StatUtil::discrete_rnd(uniformDistribution)];
                
				particleCandidates[iParticle].fromParticle = iParticle;
                particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;

                particleCandidates[iCandidate].iBestChannelOrder = processedParticle->iMaxChannelOrderAPP();
                particleCandidates[iCandidate].likelihood = 1.0;
                particleCandidates[iCandidate].weight = processedParticle->getWeight();
                normConst += particleCandidates[iCandidate].weight;
            }
            iCandidate = _particleFilter->nParticles();
        }

        // a vector of size the number of generated candidates is declared...
        VectorXd weights(iCandidate);

        // ...to store their weights
        for(uint i=0;i<iCandidate;i++)
            weights(i) = particleCandidates[i].weight/normConst;

        // the candidates that are going to give rise to particles are selected
        vector<uint> indexesSelectedCandidates = _resamplingAlgorithm->obtainIndexes(_particleFilter->capacity(),weights);

        // every survivor candidate is associated with an old particle
        vector<uint> indexesParticles(indexesSelectedCandidates.size());
        for(uint i=0;i<indexesSelectedCandidates.size();i++)
        {
            indexesParticles[i] = particleCandidates[indexesSelectedCandidates[i]].fromParticle;
            _particlesBestChannelOrders[i] = particleCandidates[indexesSelectedCandidates[i]].iBestChannelOrder;
        }

        // the chosen particles are kept without modification (yet)
        _particleFilter->keepParticles(indexesParticles);

        // every surviving particle is modified according to what it says its corresponding candidate
        for(uint iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)
        {
            processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->getParticle(iParticle));

            // sampled symbols are copied into the corresponding particle
            processedParticle->setSymbolVector(iObservationToBeProcessed,particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.col(_maxOrder-1));

            for(uint iChannelOrder=0;iChannelOrder<processedParticle->nChannelMatrixEstimators();iChannelOrder++)
            {
                // channel matrix is estimated by means of the particle channel estimator
                processedParticle->setChannelMatrix(iChannelOrder,iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(iChannelOrder)->nextMatrix(observations.col(iObservationToBeProcessed),particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.block(0,_maxOrder-_candidateOrders[iChannelOrder],_nInputs,_candidateOrders[iChannelOrder]),noiseVariances[iObservationToBeProcessed]));

				// notice that the likelihood of a candidate is the sum of the likelihoods for the different channel orders, i.e, the normalization constant
                processedParticle->setChannelOrderAPP(particleCandidates[indexesSelectedCandidates[iParticle]].unnormalizedChannelOrderAPPs(iChannelOrder)/particleCandidates[indexesSelectedCandidates[iParticle]].likelihood,iChannelOrder);
            }

            processedParticle->setWeight(particleCandidates[indexesSelectedCandidates[iParticle]].weight);

        } // for(int iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)

        _particleFilter->normalizeWeights();

        processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->getBestParticle());

        for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
            _channelOrderAPPs[0](iChannelOrder,iObservationToBeProcessed) = processedParticle->getChannelOrderAPP(iChannelOrder);

		// this doesn't make much sense...probably some debug stuff
        if(_particlesBestChannelOrders[_particleFilter->iBestParticle()]==iBestChannelOrder)
            timesBestChannelOrder++;
        else
        {
            iBestChannelOrder = _particlesBestChannelOrders[_particleFilter->iBestParticle()];
            timesBestChannelOrder = 0;
        }

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected+_d;iObservationToBeProcessed++)

    delete[] particleCandidates;
}

uint MLSDmAlgorithm::iBestChannelOrder(uint iBestParticle)
{
    return _particlesBestChannelOrders[iBestParticle];
}

void MLSDmAlgorithm::beforeInitializingParticles(const MatrixXd &observations,vector<double> &noiseVariances,const MatrixXd &trainingSequence)
{
    // in this case beforeInitializingParticles computes the APP probabilities that are obtained after the training sequence

    MatrixXd sequenceToProcess(trainingSequence.rows(),_preamble.cols()+trainingSequence.cols());
    sequenceToProcess << _preamble,trainingSequence;

    if(observations.cols() < (_iFirstObservation+trainingSequence.cols()))
        throw RuntimeException("MLSDmAlgorithm::beforeInitializingParticles: not enough number of observations.");

    vector<ChannelMatrixEstimator *> clonedChannelEstimators(_candidateOrders.size());
    for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
    {
        // the initial channel estimators are kept for computing the APP of the channel orders
        clonedChannelEstimators[iOrder] = _channelEstimators[iOrder]->clone();

        // at the beginning, all the channel orders have the same probability
        _channelOrderAPPs[0](iOrder,_preamble.cols()-1) = 1.0/double(_candidateOrders.size());
    }

    double normConst;
    vector<double> unnormalizedChannelOrderAPPs(_candidateOrders.size());

    for(int i=_preamble.cols();i<sequenceToProcess.cols();i++)
    {
        normConst = 0.0;
        for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
        {
            // unnormalized channel order APP
            unnormalizedChannelOrderAPPs[iOrder] = _channelOrderAPPs[0](iOrder,i-1)*clonedChannelEstimators[iOrder]->likelihood(observations.col(i),sequenceToProcess.block(0,i-_candidateOrders[iOrder]+1,_nInputs,_candidateOrders[iOrder]),noiseVariances[i]);
            normConst += unnormalizedChannelOrderAPPs[iOrder];

            clonedChannelEstimators[iOrder]->nextMatrix(observations.col(i),sequenceToProcess.block(0,i-_candidateOrders[iOrder]+1,_nInputs,_candidateOrders[iOrder]),noiseVariances[i]);
        }

        if(normConst!=0.0)
            for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
                _channelOrderAPPs[0](iOrder,i) = unnormalizedChannelOrderAPPs[iOrder] / normConst;
        else
            for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
                _channelOrderAPPs[0](iOrder,i) = _channelOrderAPPs[0](iOrder,i-1);
    }

    for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
        delete clonedChannelEstimators[iOrder];
}
