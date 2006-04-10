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

SMCAlgorithm::SMCAlgorithm(string name, Alphabet alphabet, ChannelMatrixEstimator& channelEstimator, tMatrix preamble,int smoothingLag,int nParticles,ResamplingCriterion resamplingCriterion): KnownChannelOrderAlgorithm(name, alphabet, channelEstimator, preamble),
// _variables initialization
_d(smoothingLag),_nParticles(nParticles),_resamplingCriterion(resamplingCriterion),_estimatedChannelMatrices(new tMatrix*[_nParticles]),_detectedSymbols(new tMatrix[_nParticles]),_particlesChannelMatrixEstimators(new ChannelMatrixEstimator*[_nParticles]),_weights(_nParticles),_reservedMemory(false)
{
// 	_estimatedChannelMatrices = new tMatrix*[_nParticles];

	// at first, we assume that all observations from the preamble need to processed
	_startDetectionTime = _m - 1;

	for(int i=0;i<_nParticles;i++)
		_particlesChannelMatrixEstimators[i] = _channelEstimator.Clone();

	_weights = 1.0/_nParticles;
}

void SMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
	int nObservations = observations.cols();
	_endDetectionTime = nObservations - _d;

	if(nObservations<(_startDetectionTime+1+_d))
		throw RuntimeException("Not enough observations.");
	
	// memory is reserved
	for(int iParticle=0;iParticle<_nParticles;iParticle++)
	{
		_estimatedChannelMatrices[iParticle] = new tMatrix[_endDetectionTime];
		_detectedSymbols[iParticle] = *(new tMatrix(_N,_endDetectionTime));
	}

	this->Process(observations,noiseVariances);
}

void SMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
	if(observations.rows()!=_L || trainingSequence.rows()!=_N)
		throw RuntimeException("Observations matrix or training sequence dimensions are wrong.");

	int iParticle,j,nObservations;

	nObservations = observations.cols();
	_endDetectionTime = nObservations - _d;

	// to process the training sequence, we need both the preamble and the symbol vectors related to it
	tMatrix preambleTrainingSequence = Util::Append(_preamble,trainingSequence);

	vector<tMatrix> trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

	for(iParticle=0;iParticle<_nParticles;iParticle++)
	{
		// memory is reserved for the current particle
		_estimatedChannelMatrices[iParticle] = new tMatrix[_endDetectionTime];
		_detectedSymbols[iParticle] = *(new tMatrix(_N,_endDetectionTime));

		//the channel estimation given by the training sequence is copied into each particle...
		for(j=_m-1;j<trainingSequenceChannelMatrices.size();j++)
		{
			_estimatedChannelMatrices[iParticle][j] = trainingSequenceChannelMatrices[j];
		}

		//... the symbols are considered detected...
		_detectedSymbols[iParticle](*(new tRange(0,_N-1)),*(new tRange(0,preambleTrainingSequence.cols()-1))).inject(preambleTrainingSequence);

		// ... and the channel estimators of all the particles are updated
		_particlesChannelMatrixEstimators[iParticle] = _channelEstimator.Clone();
	}

	// the Process method must start in
	_startDetectionTime = trainingSequenceChannelMatrices.size();

	this->Process(observations,noiseVariances);
}
