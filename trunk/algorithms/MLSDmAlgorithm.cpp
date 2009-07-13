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

#define DEBUG4

MLSDmAlgorithm::MLSDmAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance): MultipleChannelEstimatorsPerParticleSMCAlgorithm (name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_particleFilter(new ParticleFilter(nParticles)),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_particlesBestChannelOrders(nParticles)
{
}


MLSDmAlgorithm::~MLSDmAlgorithm()
{
    delete _particleFilter;
}

void MLSDmAlgorithm::InitializeParticles()
{
    vector<ChannelMatrixEstimator *> channelEstimatorsClone(_channelEstimators.size());
    for(uint i=0;i<_candidateOrders.size();i++)
        channelEstimatorsClone[i] = _channelEstimators[i]->clone();

    // we begin with only one particle
    ParticleWithChannelEstimationAndChannelOrderAPP *particle = new ParticleWithChannelEstimationAndChannelOrderAPP(1.0,_nInputs,_iLastSymbolVectorToBeDetected+_d,channelEstimatorsClone);

    particle->setSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);

    // the available APP's just before the _startDetectionTime instant are copied into the particle
    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        particle->setChannelOrderAPP(_channelOrderAPPs(iChannelOrder,_startDetectionTime-1),iChannelOrder);

    _particleFilter->AddParticle(particle);
}

void MLSDmAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
    uint nSymbolVectors = (int) pow((double)_alphabet.length(),(double)_nInputs);
    tRange rMaxChannelOrderMinus1FirstColumns(0,_maxOrder-2),rAll;
    vector<tSymbol> testedVector(_nInputs);
    tVector computedObservations(_nOutputs);
    int iCandidate,m,iBestUnnormalizedChannelOrderAPP,k,iParticle;
    uint iChannelOrder,iTestedVector;
    ParticleWithChannelEstimationAndChannelOrderAPP *processedParticle;
    double normConst,likelihood;

    typedef struct{
        int fromParticle;
        tMatrix symbolVectorsMatrix;
        int iBestChannelOrder;
        tVector unnormalizedChannelOrderAPPs;
        double likelihood;
        double weight;
    }tParticleCandidate;

    tParticleCandidate *particleCandidates = new tParticleCandidate[_particleFilter->Capacity()*nSymbolVectors];

    // "symbolVectorsMatrix" will contain all the symbols involved in the current observation
    tMatrix symbolVectorsMatrix(_nInputs,_maxOrder);
    tVector symbolsVector;

    int lastSymbolVectorStart = _nInputsXchannelOrderaxOrder - _nInputs;

    tRange rMaxChannelOrderMinus1PrecedentColumns(_startDetectionTime-_maxOrder+1,_startDetectionTime-1);

    vector<bool> activeCandidateOrders(_candidateOrders.size(),true);
    int iBestChannelOrder = 0,timesBestChannelOrder = 0;

    for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected+_d;iObservationToBeProcessed++)
    {
#ifdef DEBUG
//      cout << iObservationToBeProcessed << endl;
#endif
        // it keeps track of the place where a new tParticleCandidate will be stored within the array
        iCandidate = 0;

        normConst = 0.0;

        // the candidates from all the particles are generated
        for(iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
        {
            processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->GetParticle(iParticle));

            symbolVectorsMatrix(rAll,rMaxChannelOrderMinus1FirstColumns).inject(processedParticle->getSymbolVectors(rMaxChannelOrderMinus1PrecedentColumns));
            symbolsVector = Util::toVector(symbolVectorsMatrix,columnwise);

            for(iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
            {
                // the corresponding testing vector is generated from the index
                _alphabet.int2symbolsArray(iTestedVector,testedVector);

                // current tested vector is copied in the m-th position
                for(k=0;k<_nInputs;k++)
                    symbolVectorsMatrix(k,_maxOrder-1) = symbolsVector(lastSymbolVectorStart+k) = testedVector[k];

                likelihood = 0.0;

                particleCandidates[iCandidate].unnormalizedChannelOrderAPPs = tVector(_candidateOrders.size());

                for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
                {
                    m = _candidateOrders[iChannelOrder];

                    tMatrix involvedSymbolVectors = symbolVectorsMatrix(rAll,tRange(_maxOrder-m,_maxOrder-1)).copy();

                    // the AR coefficiente is accounted for
//                  involvedSymbolVectors *= _ARcoefficient; <--------------------------------- (when Kalman estimator, it is already accounted for)

                    particleCandidates[iCandidate].unnormalizedChannelOrderAPPs(iChannelOrder) = processedParticle->getChannelOrderAPP(iChannelOrder)*processedParticle->getChannelMatrixEstimator(iChannelOrder)->likelihood(observations.col(iObservationToBeProcessed),involvedSymbolVectors,noiseVariances[iObservationToBeProcessed]);


                    likelihood += particleCandidates[iCandidate].unnormalizedChannelOrderAPPs(iChannelOrder);
                }

                // if the channelOrderNormConst is zero, we don't generate a candidate for this particle and this symbol vector
                if(likelihood==0.0)
                    continue;
//                  throw RuntimeException("UTSAlgorithm::Process: likelihood is zero.");

                Util::max(particleCandidates[iCandidate].unnormalizedChannelOrderAPPs,iBestUnnormalizedChannelOrderAPP);

                particleCandidates[iCandidate].fromParticle = iParticle;
                particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;
                particleCandidates[iCandidate].iBestChannelOrder = iBestUnnormalizedChannelOrderAPP;
                particleCandidates[iCandidate].likelihood = likelihood;
                particleCandidates[iCandidate].weight = processedParticle->getWeight()*likelihood;
                normConst += particleCandidates[iCandidate].weight;

                iCandidate++;
            } // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

        } // for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)

        // if none of the candidates was valid
        if(iCandidate==0)
        {
            tVector uniformDistribution(_alphabet.length());
            for(int iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                uniformDistribution(iAlphabet) = 1.0/_alphabet.length();

            for(iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
            {
                processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->GetParticle(iParticle));
                symbolVectorsMatrix(rAll,rMaxChannelOrderMinus1FirstColumns).inject(processedParticle->getSymbolVectors(rMaxChannelOrderMinus1PrecedentColumns));
                for(k=0;k<_nInputs;k++)
                    symbolVectorsMatrix(k,_maxOrder-1) = _alphabet[StatUtil::discrete_rnd(uniformDistribution)];
                particleCandidates[iParticle].fromParticle = iParticle;
                particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;

                particleCandidates[iCandidate].iBestChannelOrder = processedParticle->iMaxChannelOrderAPP();
                particleCandidates[iCandidate].likelihood = 1.0;
                particleCandidates[iCandidate].weight = processedParticle->getWeight();
                normConst += particleCandidates[iCandidate].weight;
            }
            iCandidate = _particleFilter->Nparticles();
        }

        // a vector of size the number of generated candidates is declared...
        tVector weights(iCandidate);

        // ...to store their weights
        for(int i=0;i<iCandidate;i++)
            weights(i) = particleCandidates[i].weight/normConst;

        // the candidates that are going to give particles are selected
        vector<int> indexesSelectedCandidates = _resamplingAlgorithm->ObtainIndexes(_particleFilter->Capacity(),weights);

        // every survivor candidate is associated with an old particle
        vector<int> indexesParticles(indexesSelectedCandidates.size());
        for(uint i=0;i<indexesSelectedCandidates.size();i++)
        {
            indexesParticles[i] = particleCandidates[indexesSelectedCandidates[i]].fromParticle;
            _particlesBestChannelOrders[i] = particleCandidates[indexesSelectedCandidates[i]].iBestChannelOrder;
        }

        // the chosen particles are kept without modification (yet)
        _particleFilter->keepParticles(indexesParticles);

        // every surviving particle is modified according to what it says its corresponding candidate
        for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
        {
            processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->GetParticle(iParticle));

            // sampled symbols are copied into the corresponding particle
            processedParticle->setSymbolVector(iObservationToBeProcessed,particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.col(_maxOrder-1));

            for(int iChannelOrder=0;iChannelOrder<processedParticle->nChannelMatrixEstimators();iChannelOrder++)
            {
                // channel matrix is estimated by means of the particle channel estimator
                processedParticle->setChannelMatrix(iChannelOrder,iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(iChannelOrder)->nextMatrix(observations.col(iObservationToBeProcessed),particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix(rAll,tRange(_maxOrder-_candidateOrders[iChannelOrder],_maxOrder-1)),noiseVariances[iObservationToBeProcessed]));

                processedParticle->setChannelOrderAPP(particleCandidates[indexesSelectedCandidates[iParticle]].unnormalizedChannelOrderAPPs(iChannelOrder)/particleCandidates[indexesSelectedCandidates[iParticle]].likelihood,iChannelOrder);
            }

            processedParticle->setWeight(particleCandidates[indexesSelectedCandidates[iParticle]].weight);

        } // for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)

        _particleFilter->NormalizeWeights();

        processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrderAPP *> (_particleFilter->GetBestParticle());

        for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
            _channelOrderAPPs(iChannelOrder,iObservationToBeProcessed) = processedParticle->getChannelOrderAPP(iChannelOrder);

        if(_particlesBestChannelOrders[_particleFilter->iBestParticle()]==iBestChannelOrder)
            timesBestChannelOrder++;
        else
        {
            iBestChannelOrder = _particlesBestChannelOrders[_particleFilter->iBestParticle()];
            timesBestChannelOrder = 0;
        }

        rMaxChannelOrderMinus1PrecedentColumns = rMaxChannelOrderMinus1PrecedentColumns + 1;

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected+_d;iObservationToBeProcessed++)

    delete[] particleCandidates;
}

int MLSDmAlgorithm::BestChannelOrderIndex(int iBestParticle)
{
    return _particlesBestChannelOrders[iBestParticle];
}

// vector<vector<tMatrix> > MLSDmAlgorithm::EstimateChannelFromTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
// {
//     tMatrix sequenceToProcess = Util::append(_preamble,trainingSequence);
//
//     if(observations.cols() < (_iFirstObservation+trainingSequence.cols()))
//         throw RuntimeException("MLSDmAlgorithm::EstimateChannelFromTrainingSequence: Insufficient number of observations.");
//
//     vector<ChannelMatrixEstimator *> initialChannelEstimators(_candidateOrders.size());
//     for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
//     {
//         // the initial channel estimators are kept for computing the APP of the channel orders
//         initialChannelEstimators[iOrder] = _channelEstimators[iOrder]->clone();
//
//         // at the beginning, all the channel orders have the same probability
//         _channelOrderAPPs(iOrder,_preamble.cols()-1) = 1.0/double(_candidateOrders.size());
//     }
//
//     tRange rAll;
//     double normConst;
//     vector<double> unnormalizedChannelOrderAPPs(_candidateOrders.size());
//
//     for(int i=_preamble.cols();i<sequenceToProcess.cols();i++)
//     {
//         normConst = 0.0;
//         for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
//         {
//             // unnormalized channel order APP
//             unnormalizedChannelOrderAPPs[iOrder] = _channelOrderAPPs(iOrder,i-1)*initialChannelEstimators[iOrder]->likelihood(observations.col(i),sequenceToProcess(rAll,tRange(i-_candidateOrders[iOrder]+1,i)),noiseVariances[i]);
//             normConst += unnormalizedChannelOrderAPPs[iOrder];
//
//             initialChannelEstimators[iOrder]->nextMatrix(observations.col(i),sequenceToProcess(rAll,tRange(i-_candidateOrders[iOrder]+1,i)),noiseVariances[i]);
//         }
//
//         if(normConst!=0.0)
//             for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
//                 _channelOrderAPPs(iOrder,i) = unnormalizedChannelOrderAPPs[iOrder] / normConst;
//         else
//             for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
//                 _channelOrderAPPs(iOrder,i) = _channelOrderAPPs(iOrder,i-1);
//     }
//
//     for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
//         delete initialChannelEstimators[iOrder];
//
//     return UnknownChannelOrderAlgorithm::EstimateChannelFromTrainingSequence(observations,noiseVariances,trainingSequence);
// }

void MLSDmAlgorithm::BeforeInitializingParticles(const tMatrix &observations,vector<double> &noiseVariances,const tMatrix &trainingSequence)
{
    // in this case BeforeInitializingParticles computes the APP probabilities that are obtained after the training sequence

    tMatrix sequenceToProcess = Util::append(_preamble,trainingSequence);

    if(observations.cols() < (_iFirstObservation+trainingSequence.cols()))
        throw RuntimeException("BeforeInitializingParticles::BeforeInitializingParticles: Insufficient number of observations.");

    vector<ChannelMatrixEstimator *> clonedChannelEstimators(_candidateOrders.size());
    for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
    {
        // the initial channel estimators are kept for computing the APP of the channel orders
        clonedChannelEstimators[iOrder] = _channelEstimators[iOrder]->clone();

        // at the beginning, all the channel orders have the same probability
        _channelOrderAPPs(iOrder,_preamble.cols()-1) = 1.0/double(_candidateOrders.size());
    }

    tRange rAll;
    double normConst;
    vector<double> unnormalizedChannelOrderAPPs(_candidateOrders.size());

    for(int i=_preamble.cols();i<sequenceToProcess.cols();i++)
    {
        normConst = 0.0;
        for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
        {
            // unnormalized channel order APP
            unnormalizedChannelOrderAPPs[iOrder] = _channelOrderAPPs(iOrder,i-1)*clonedChannelEstimators[iOrder]->likelihood(observations.col(i),sequenceToProcess(rAll,tRange(i-_candidateOrders[iOrder]+1,i)),noiseVariances[i]);
            normConst += unnormalizedChannelOrderAPPs[iOrder];

            clonedChannelEstimators[iOrder]->nextMatrix(observations.col(i),sequenceToProcess(rAll,tRange(i-_candidateOrders[iOrder]+1,i)),noiseVariances[i]);
        }

        if(normConst!=0.0)
            for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
                _channelOrderAPPs(iOrder,i) = unnormalizedChannelOrderAPPs[iOrder] / normConst;
        else
            for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
                _channelOrderAPPs(iOrder,i) = _channelOrderAPPs(iOrder,i-1);
    }

    for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
        delete clonedChannelEstimators[iOrder];
}
