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
#include "MultipleChannelEstimatorsPerParticleSMCAlgorithm.h"

MultipleChannelEstimatorsPerParticleSMCAlgorithm::MultipleChannelEstimatorsPerParticleSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation,int smoothingLag,int nParticles,ResamplingAlgorithm *resamplingAlgorithm): UnknownChannelOrderAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation)
//variables initialization
,_d(smoothingLag),_allSymbolsRows(0,_N-1),_resamplingAlgorithm(resamplingAlgorithm),_particleFilter(nParticles,_candidateOrders),_nParticlesPerChannelOrder(_candidateOrders.size())
{
    // at first, we assume that all symbol vectors from the preamble need to be processed
    _startDetectionTime = _preamble.cols();
}


MultipleChannelEstimatorsPerParticleSMCAlgorithm::~MultipleChannelEstimatorsPerParticleSMCAlgorithm()
{
}

void MultipleChannelEstimatorsPerParticleSMCAlgorithm::InitializeParticles()
{
    int iParticlePresentOrder,nParticlesPresentOrder,iParticle=0;
    int iChannelMatrixEstimator;

    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
    {
    // If the number of particles is not divisible by the number of channel orders, the remaining particles are assigned to the first channel orders (the + (......) term)
        nParticlesPresentOrder = _particleFilter.Nparticles()/_candidateOrders.size() + (_particleFilter.Nparticles() % _candidateOrders.size() > iChannelOrder);

        for(iParticlePresentOrder=0;iParticlePresentOrder<nParticlesPresentOrder;iParticlePresentOrder++,iParticle++)
        {
            // a clone of each of the channel matrix estimators is constructed...
            vector< ChannelMatrixEstimator * > thisParticleChannelMatrixEstimators(_candidateOrders.size());
            for(iChannelMatrixEstimator=0;iChannelMatrixEstimator<_candidateOrders.size();iChannelMatrixEstimator++)
                thisParticleChannelMatrixEstimators[iChannelMatrixEstimator] = _channelEstimators[iChannelMatrixEstimator]->Clone();

            // ... and passed within a vector to each particle
            _particleFilter.SetParticle(new ParticleWithChannelEstimationAndChannelOrder(1.0/(double)_particleFilter.Nparticles(),_N,_K,thisParticleChannelMatrixEstimators,_candidateOrders[iChannelOrder]),iParticle);
        }
    }
}

void MultipleChannelEstimatorsPerParticleSMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();

    if(nObservations<_startDetectionTime+_maxOrder)
        throw RuntimeException("MultipleChannelEstimatorsPerParticleSMCAlgorithm::Run: Not enough observations.");

    this->InitializeParticles();

    this->Process(observations,noiseVariances);
}

void MultipleChannelEstimatorsPerParticleSMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
    if(observations.rows()!=_L || trainingSequence.rows()!=_N)
        throw RuntimeException("MultipleChannelEstimatorsPerParticleSMCAlgorithm::Run: Observations matrix or training sequence dimensions are wrong.");

    // observations are going to be needed to find the best particle
//     _observations = observations;

//     int iParticle,iParticlePresentOrder,j;
    int iParticle,j,iChannelOrder;

    // to process the training sequence, we need both the preamble and the symbol vectors related to it
    tMatrix preambleTrainingSequence = Util::Append(_preamble,trainingSequence);


    tRange rSymbolVectorsTrainingSequece(0,preambleTrainingSequence.cols()-1);

    vector<vector<tMatrix> > trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

    this->InitializeParticles();

//     for(int iChannelOrder=0,iParticle=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
//     {
//         for(iParticlePresentOrder=0;iParticlePresentOrder<_nParticlesPerChannelOrder[iChannelOrder];iParticlePresentOrder++,iParticle++)
//         {
//             ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particleFilter.GetParticle(iParticle));
//
//             //the channel estimation given by the training sequence is copied into each particle...
//             for(j=0;j<trainingSequenceChannelMatrices[iChannelOrder].size();j++)
//             {
//                 processedParticle->SetChannelMatrix(_preamble.cols()+j,trainingSequenceChannelMatrices[iChannelOrder][j]);
//             }
//
//             //... the symbols are considered detected...
//             processedParticle->SetSymbolVectors(rSymbolVectorsTrainingSequece,preambleTrainingSequence);
//         }
//     }

    for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
    {
        ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particleFilter.GetParticle(iParticle));

        for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        {
            //the channel estimation given by the training sequence is copied into each particle...
            for(j=0;j<trainingSequenceChannelMatrices[iChannelOrder].size();j++)
            {
                processedParticle->SetChannelMatrix(iChannelOrder,_preamble.cols()+j,trainingSequenceChannelMatrices[iChannelOrder][j]);
            }
        }

        //... the symbols are considered detected...
        processedParticle->SetSymbolVectors(rSymbolVectorsTrainingSequece,preambleTrainingSequence);
    }

    // the Process method must start in
    _startDetectionTime += trainingSequence.cols();

    this->Process(observations,noiseVariances);
}

tMatrix MultipleChannelEstimatorsPerParticleSMCAlgorithm::GetDetectedSymbolVectors()
{
    int iBestParticle;
    Util::Max(_particleFilter.GetWeightsVector(),iBestParticle);

    ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrder *> ( _particleFilter.GetParticle(iBestParticle));

    cout << "La particula seleccionada tiene orden " << processedParticle->GetChannelOrder() << endl;

    return ((_particleFilter.GetParticle(iBestParticle))->GetAllSymbolVectors())(_allSymbolsRows,tRange(_preamble.cols(),_K-1));
}

vector<tMatrix> MultipleChannelEstimatorsPerParticleSMCAlgorithm::GetEstimatedChannelMatrices()
{
    vector<tMatrix> channelMatrices(0);
    return channelMatrices;
}
