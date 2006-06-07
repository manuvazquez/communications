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
#include "SMCAlgorithm.h"

SMCAlgorithm::SMCAlgorithm(string name, Alphabet alphabet, ChannelMatrixEstimator *channelEstimator, tMatrix preamble,int smoothingLag,int nParticles,ResamplingCriterion resamplingCriterion,StdResamplingAlgorithm resamplingAlgorithm): KnownChannelOrderAlgorithm(name, alphabet, channelEstimator, preamble),
// _variables initialization
_d(smoothingLag),_nParticles(nParticles),_resamplingCriterion(resamplingCriterion),_resamplingAlgorithm(resamplingAlgorithm),
_particles(new ParticleWithChannelEstimation*[nParticles]),
_allSymbolsRows(0,_N-1)
{
	// at first, we assume that all observations from the preamble need to be processed
	_startDetectionTime = _m - 1;

	for(int i=0;i<_nParticles;i++)
	{
		_particles[i] = NULL;
	}
}

SMCAlgorithm::~SMCAlgorithm()
{
	for(int i=0;i<_nParticles;i++)
	{
		delete _particles[i];
	}

	delete[] _particles;
}

void SMCAlgorithm::InitializeParticles()
{
	// memory is reserved
	for(int iParticle=0;iParticle<_nParticles;iParticle++)
	{
		_particles[iParticle] = new ParticleWithChannelEstimation(1.0/(double)_nParticles,_N,_endDetectionTime,_channelEstimator->Clone());
	}
}

void SMCAlgorithm::Run(const tMatrix &observations,vector<double> noiseVariances)
{
	int nObservations = observations.cols();
	_endDetectionTime = nObservations - _d;

	if(nObservations<(_startDetectionTime+1+_d))
		throw RuntimeException("SMCAlgorithm::Run: Not enough observations.");

	this->InitializeParticles();

	this->Process(observations,noiseVariances);
}

void SMCAlgorithm::Run(const tMatrix &observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
// 	cout << "Running with training sequence..." << endl;

	if(observations.rows()!=_L || trainingSequence.rows()!=_N)
		throw RuntimeException("Run: Observations matrix or training sequence dimensions are wrong.");

	int iParticle,j,nObservations;

	nObservations = observations.cols();
	_endDetectionTime = nObservations - _d;

	// to process the training sequence, we need both the preamble and the symbol vectors related to it
	tMatrix preambleTrainingSequence = Util::Append(_preamble,trainingSequence);


	tRange rSymbolVectorsTrainingSequece(0,preambleTrainingSequence.cols()-1);

	vector<tMatrix> trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

	this->InitializeParticles();

	for(iParticle=0;iParticle<_nParticles;iParticle++)
	{
		//the channel estimation given by the training sequence is copied into each particle...
		for(j=_m-1;j<trainingSequenceChannelMatrices.size();j++)
		{
// 			_estimatedChannelMatrices[iParticle][j] = trainingSequenceChannelMatrices[j];
			_particles[iParticle]->SetChannelMatrix(j,trainingSequenceChannelMatrices[j]);
		}

		//... the symbols are considered detected...
		_particles[iParticle]->SetSymbolVectors(rSymbolVectorsTrainingSequece,preambleTrainingSequence);
	}

	// the Process method must start in
	_startDetectionTime = trainingSequenceChannelMatrices.size();

	this->Process(observations,noiseVariances);
}

void SMCAlgorithm::Resampling(int endResamplingTime)
{
	if(_resamplingCriterion.ResamplingNeeded(_particles,_nParticles))
	{
		vector<int> indexes = StatUtil::Discrete_rnd(_nParticles,GetWeightsVector());

		_resamplingAlgorithm.Resampling(&_particles,_nParticles,indexes);

// 		cout << "Resampling..." << endl;
	}
}

tMatrix SMCAlgorithm::GetDetectedSymbolVectors()
{
    // best particle is chosen
    int iBestParticle;
    Util::Max(GetWeightsVector(),iBestParticle);

    return (_particles[iBestParticle]->GetAllSymbolVectors())(_allSymbolsRows,tRange(_m-1,_endDetectionTime-1));
}

vector<tMatrix> SMCAlgorithm::GetDetectedChannelMatrices()
{
    vector<tMatrix> channelMatrices;
    channelMatrices.reserve(_endDetectionTime-_m+1);

    // best particle is chosen
    int iBestParticle;
    Util::Max(GetWeightsVector(),iBestParticle);
    
    for(int i=_m-1;i<_endDetectionTime;i++)
        channelMatrices.push_back(_particles[iBestParticle]->GetChannelMatrix(i));

    return channelMatrices;
}
