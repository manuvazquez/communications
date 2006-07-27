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
	// at first, we assume that all symbol vectors from the preamble need to be processed
    _startDetectionSymbolVector = _preamble.cols();

	// and all the observations from _firstObservationIndex
	_startDetectionObservation = _firstObservationIndex;

	// we set the number of particles that will be assignated to each channel order
	for(int iChannelOrder=0;iChannelOrder<_nCandidateOrders;iChannelOrder++)
	{
		// the particles are distributed among the different channel orders. If the number of particles is not divisible by the number of channel orders, the remaining particles are assigned to the first channel orders (the + (......) term)
		_nParticlesPerChannelOrder[iChannelOrder] = _particleFilter.Nparticles()/_nCandidateOrders + (_particleFilter.Nparticles() % _nCandidateOrders > iChannelOrder);
	}

	_channelOrderWeightsSum = new double[_maxOrder+1];
}


UnknownChannelOrderSMCAlgorithm::~UnknownChannelOrderSMCAlgorithm()
{
	delete[] _channelOrderWeightsSum;
}

void UnknownChannelOrderSMCAlgorithm::InitializeParticles()
{
	// we have to initialize _nParticlesPerChannelOrder[i] particles of order _candidateOrders[i]
	int iParticlePresentOrder,iParticle=0;
	for(int iChannelOrder=0;iChannelOrder<_nCandidateOrders;iChannelOrder++)
		for(iParticlePresentOrder=0;iParticlePresentOrder<_nParticlesPerChannelOrder[iChannelOrder];iParticlePresentOrder++,iParticle++)
		{
        	_particleFilter.SetParticle(new ParticleWithChannelEstimationAndChannelOrder(1.0/(double)_particleFilter.Nparticles(),_N,_K-_firstObservationIndex+_preamble.cols(),_channelEstimators[iChannelOrder]->Clone(),_candidateOrders[iChannelOrder]),iParticle);
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
    if(observations.rows()!=_L || trainingSequence.rows()!=_N)
        throw RuntimeException("Run: Observations matrix or training sequence dimensions are wrong.");

	// observations are going to be needed to find the best particle
	_observations = observations;

    int iParticle,iParticlePresentOrder,j;

    // to process the training sequence, we need both the preamble and the symbol vectors related to it
    tMatrix preambleTrainingSequence = Util::Append(_preamble,trainingSequence);


    tRange rSymbolVectorsTrainingSequece(0,preambleTrainingSequence.cols()-1);

    vector<vector<tMatrix> > trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

    this->InitializeParticles();

	for(int iChannelOrder=0,iParticle=0;iChannelOrder<_nCandidateOrders;iChannelOrder++)
    {
		for(iParticlePresentOrder=0;iParticlePresentOrder<_nParticlesPerChannelOrder[iChannelOrder];iParticlePresentOrder++,iParticle++)
		{
			ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particleFilter.GetParticle(iParticle));
	
			//the channel estimation given by the training sequence is copied into each particle...
			for(j=0;j<trainingSequenceChannelMatrices[iChannelOrder].size();j++)
			{
				processedParticle->SetChannelMatrix(_preamble.cols()+j,trainingSequenceChannelMatrices[iChannelOrder][j]);
			}
	
			//... the symbols are considered detected...
			processedParticle->SetSymbolVectors(rSymbolVectorsTrainingSequece,preambleTrainingSequence);
		}
    }

    // the Process method must start in
    _startDetectionObservation = _firstObservationIndex + trainingSequence.cols();
	_startDetectionSymbolVector = _preamble.cols() + trainingSequence.cols();

    this->Process(observations,noiseVariances);
}

void UnknownChannelOrderSMCAlgorithm::NormalizeParticleGroups()
{
	int iParticle;

	// the sum of the weight of each group of particles is set to zero...
	for(int iChannelOrder=0;iChannelOrder<_nCandidateOrders;iChannelOrder++)
		_channelOrderWeightsSum[_candidateOrders[iChannelOrder]] = 0.0;
	
	// ... to compute it here
	for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
	{
		ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particleFilter.GetParticle(iParticle));

		_channelOrderWeightsSum[processedParticle->GetChannelOrder()] += processedParticle->GetWeight();
	}

	// each particle is normalized according to its channel order
	for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
	{
		ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particleFilter.GetParticle(iParticle));

		processedParticle->SetWeight(processedParticle->GetWeight()/_channelOrderWeightsSum[processedParticle->GetChannelOrder()]);
	}	
}

void UnknownChannelOrderSMCAlgorithm::Resampling()
{
// 	cout << "Antes de resampling quedan:" << endl;
// 	vector<vector<int> > indices2 = GetIndexesOfChannelOrders();
// 	for(int i=0;i<indices2.size();i++)
// 		cout << indices2[i].size() << "particulas de orden " << i << endl;

	tVector weigths = _particleFilter.GetWeightsVector();

    if(_particleFilter._resamplingCriterion.ResamplingNeeded(weigths))
    {
        vector<int> indexes = StatUtil::Discrete_rnd(_particleFilter.Nparticles(),weigths);
		_particleFilter.SelectParticles(indexes);
    }
}

void UnknownChannelOrderSMCAlgorithm::ResamplingByParticleGroups()
{
	vector<vector <int> > indexesOfChannelOrders = GetIndexesOfChannelOrders();
	tVector weights = _particleFilter.GetWeightsVector();

	for(int iChannelOrder=0;iChannelOrder<_nCandidateOrders;iChannelOrder++)
	{
		// if there is no particles with this channel order
		if(_nParticlesPerChannelOrder[iChannelOrder]==0)
			continue;

		if(_particleFilter._resamplingCriterion.ResamplingNeeded(weights,indexesOfChannelOrders[iChannelOrder]))
		{
			// the weights corresponding to the channel order being processed are selected. Note that they all are adjacent
			tRange rWeightsOrder(indexesOfChannelOrders[iChannelOrder][0],indexesOfChannelOrders[iChannelOrder][indexesOfChannelOrders[iChannelOrder].size()-1]);

// 			cout << rWeightsOrder << endl;

			vector<int> auxIndexes = StatUtil::Discrete_rnd(_nParticlesPerChannelOrder[iChannelOrder],weights(rWeightsOrder));

			vector<int> indexes(_nParticlesPerChannelOrder[iChannelOrder]);
			for(int i=0;i<auxIndexes.size();i++)
					indexes[i] = indexesOfChannelOrders[iChannelOrder][auxIndexes[i]];

			_particleFilter.SelectParticles(indexes,indexesOfChannelOrders[iChannelOrder]);
		}	
	}
}

/**
 * It returns a vector of vectors of int, such that vector[i] contains the indexes of the particles of order _candidateOrders[i]
 * @return 
 */
vector<vector<int> > UnknownChannelOrderSMCAlgorithm::GetIndexesOfChannelOrders()
{
	vector<vector<int> > res(_nCandidateOrders);

	for(int iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
	{
		ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndChannelOrder *> (_particleFilter.GetParticle(iParticle));

		res[_channelOrder2index[processedParticle->GetChannelOrder()]].push_back(iParticle);
	}
	return res;
}

int UnknownChannelOrderSMCAlgorithm::BestParticle()
{
	int m,res;
	tVector computedObservation(_L), realObservation(_L);
	tVector errorVector(_L);
	double norm;
	int *iMaxWeights = new int[_nCandidateOrders];
// 	tVector particleError(_particleFilter.Nparticles());
	tVector particleError(_nCandidateOrders);
	int iMinParticleError;

	vector<vector<int> > indexesOfChannelOrders = GetIndexesOfChannelOrders();
	tVector weights = _particleFilter.GetWeightsVector();
	for(int iChannelOrder=0;iChannelOrder<_nCandidateOrders;iChannelOrder++)
	{
		// if there is no particles with this channel order
		if(_nParticlesPerChannelOrder[iChannelOrder]==0)
			continue;
		
		tRange rWeightsOrder(indexesOfChannelOrders[iChannelOrder][0],indexesOfChannelOrders[iChannelOrder][indexesOfChannelOrders[iChannelOrder].size()-1]);

		Util::Max(weights(rWeightsOrder),iMaxWeights[iChannelOrder]);

		// the best particle for this channel order is chosen
		ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrder *> ( _particleFilter.GetParticle(indexesOfChannelOrders[iChannelOrder][iMaxWeights[iChannelOrder]]));

		m = processedParticle->GetChannelOrder();

		for(int iObservationToBeProcessed=_startDetectionObservation;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
		{
			tMatrix symbols = processedParticle->GetSymbolVectors(_startDetectionSymbolVector-_startDetectionObservation+iObservationToBeProcessed-m+1,_startDetectionSymbolVector-_startDetectionObservation+iObservationToBeProcessed);

			tVector symbolsVector = Util::ToVector(symbols,columnwise);

			tMatrix channelMatrix = processedParticle->GetChannelMatrix(_startDetectionSymbolVector-_startDetectionObservation+iObservationToBeProcessed);

			// computedObservation = channelMatrix*symbols
			Blas_Mat_Vec_Mult(channelMatrix,symbolsVector,computedObservation);

			// needed because Util::Add receives a reference
			realObservation = _observations.col(iObservationToBeProcessed);

			// errorVector = realObservation - computedObservation
			Util::Add(realObservation,computedObservation,errorVector,1.0,-1.0);

			norm = Blas_Norm2(errorVector);

			particleError(iChannelOrder) += norm*norm;
		}
	}
	Util::Min(particleError,res);

	delete[] iMaxWeights;
	return res;
}

tMatrix UnknownChannelOrderSMCAlgorithm::GetDetectedSymbolVectors()
{
//     // best particle is chosen
//     int iBestParticle = BestParticle();

    int iBestParticle;
    Util::Max(_particleFilter.GetWeightsVector(),iBestParticle);

	ParticleWithChannelEstimationAndChannelOrder *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndChannelOrder *> ( _particleFilter.GetParticle(iBestParticle));

	cout << "La particula seleccionada tiene orden " << processedParticle->GetChannelOrder() << endl;

    return ((_particleFilter.GetParticle(iBestParticle))->GetAllSymbolVectors())(_allSymbolsRows,tRange(_preamble.cols(),_K-_firstObservationIndex+_preamble.cols()-1));
}

vector<tMatrix> UnknownChannelOrderSMCAlgorithm::GetEstimatedChannelMatrices()
{
//     vector<tMatrix> channelMatrices;
//     channelMatrices.reserve(_K-_preamble.cols());
// 
//     // best particle is chosen
//     int iBestParticle = BestParticle();
//     
//     for(int i=_preamble.cols();i<_K-_firstObservationIndex+_preamble.cols();i++)
//         channelMatrices.push_back((_particleFilter.GetParticle(iBestParticle))->GetChannelMatrix(i));

    vector<tMatrix> channelMatrices(0);
    return channelMatrices;
}
