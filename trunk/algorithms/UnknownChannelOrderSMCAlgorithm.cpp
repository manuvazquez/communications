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
#include "UnknownChannelOrderSMCAlgorithm.h"

UnknownChannelOrderSMCAlgorithm::UnknownChannelOrderSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int firstObservationIndex,int smoothingLag,int nParticles,ResamplingCriterion resamplingCriterion,StdResamplingAlgorithm resamplingAlgorithm): UnknownChannelOrderAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, firstObservationIndex),
//variables initialization
_d(smoothingLag),_allSymbolsRows(0,_N-1),_particleFilter(nParticles,resamplingCriterion,resamplingAlgorithm),_nParticlesPerChannelOrder(_nCandidateOrders)
{
    // at first, we assume that all observations from the preamble need to be processed
    _startDetectionTime = _preamble.cols();

	// we set the number of particles that will be assignated to each channel order
	for(int iChannelOrder=0;iChannelOrder<_nCandidateOrders;iChannelOrder++)
	{
		// the particles are distributed among the different channel orders. If the number of particles is not divisible by the number of channel orders, the remaining particles are assigned to the first channel orders (the + (......) term)
		_nParticlesPerChannelOrder[iChannelOrder] = _particleFilter.Nparticles()/_nCandidateOrders + (_particleFilter.Nparticles() % _nCandidateOrders > iChannelOrder);
	}
}


UnknownChannelOrderSMCAlgorithm::~UnknownChannelOrderSMCAlgorithm()
{
}

void UnknownChannelOrderSMCAlgorithm::InitializeParticles()
{
	// we have to initialize _nParticlesPerChannelOrder[i] particles of order _candidateOrders[i]
	int iParticlePresentOrder,iParticle=0;
	for(int iChannelOrder=0;iChannelOrder<_nCandidateOrders;iChannelOrder++)
		for(iParticlePresentOrder=0;iParticlePresentOrder<_nParticlesPerChannelOrder[iChannelOrder];iParticlePresentOrder++,iParticle++)
		{
        	_particleFilter.SetParticle(new ParticleWithChannelEstimationAndChannelOrder(1.0/(double)_particleFilter.Nparticles(),_N,_K,_channelEstimators[iChannelOrder]->Clone(),_candidateOrders[iChannelOrder]),iParticle);
		}

	for(int i=0;i<_particleFilter.Nparticles();i++)
	{
		ParticleWithChannelEstimation *part = _particleFilter.GetParticle(i);
		cout << "Particula " << i << endl << (part->GetChannelMatrixEstimator())->LastEstimatedChannelMatrix() << endl;
	}
}

void UnknownChannelOrderSMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();

    if(nObservations<_K)
        throw RuntimeException("UnknownChannelOrderSMCAlgorithm::Run: Not enough observations.");

    this->InitializeParticles();

    this->Process(observations,noiseVariances);
}

void UnknownChannelOrderSMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
//     if(observations.rows()!=_L || trainingSequence.rows()!=_N)
//         throw RuntimeException("Run: Observations matrix or training sequence dimensions are wrong.");
// 
//     int iParticle,j;
// 
//     // to process the training sequence, we need both the preamble and the symbol vectors related to it
//     tMatrix preambleTrainingSequence = Util::Append(_preamble,trainingSequence);
// 
// 
//     tRange rSymbolVectorsTrainingSequece(0,preambleTrainingSequence.cols()-1);
// 
//     vector<tMatrix> trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);
// 
//     this->InitializeParticles();
// 
//     for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
//     {
// 		ParticleWithChannelEstimation *processedParticle = _particleFilter.GetParticle(iParticle);
// 	
//         //the channel estimation given by the training sequence is copied into each particle...
//         for(j=_m-1;j<trainingSequenceChannelMatrices.size();j++)
//         {
//             processedParticle->SetChannelMatrix(j,trainingSequenceChannelMatrices[j]);
//         }
// 
//         //... the symbols are considered detected...
//         processedParticle->SetSymbolVectors(rSymbolVectorsTrainingSequece,preambleTrainingSequence);
//     }
// 
//     // the Process method must start in
//     _startDetectionTime = trainingSequenceChannelMatrices.size();
// 
//     this->Process(observations,noiseVariances);
}
