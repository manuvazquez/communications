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

SMCAlgorithm::SMCAlgorithm(string name, Alphabet alphabet, ChannelMatrixEstimator& channelEstimator, tMatrix preamble,int smoothingLag,int nParticles,ResamplingCriterion resamplingCriterion,StdResamplingAlgorithm resamplingAlgorithm): KnownChannelOrderAlgorithm(name, alphabet, channelEstimator, preamble),
// _variables initialization
_d(smoothingLag),_nParticles(nParticles),_resamplingCriterion(resamplingCriterion),_resamplingAlgorithm(resamplingAlgorithm),_estimatedChannelMatrices(new tMatrix*[_nParticles]),_detectedSymbols(new tMatrix[_nParticles]),_particlesChannelMatrixEstimators(new ChannelMatrixEstimator*[_nParticles]),_weights(_nParticles),_reservedMemory(false)
{
// 	_estimatedChannelMatrices = new tMatrix*[_nParticles];

	// at first, we assume that all observations from the preamble need to be processed
	_startDetectionTime = _m - 1;

	for(int i=0;i<_nParticles;i++)
	{
		_particlesChannelMatrixEstimators[i] = _channelEstimator.Clone();
		_estimatedChannelMatrices[i] = NULL;
	}

	_weights = 1.0/(double)_nParticles;
}

SMCAlgorithm::~ SMCAlgorithm()
{
	for(int i=0;i<_nParticles;i++)
	{
		delete _particlesChannelMatrixEstimators[i];
		delete[] _estimatedChannelMatrices[i]; // <----------------
	}
	delete[] _detectedSymbols;
	delete[] _particlesChannelMatrixEstimators;
	delete[] _estimatedChannelMatrices;
}

void SMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
	int nObservations = observations.cols();
	_endDetectionTime = nObservations - _d;

	if(nObservations<(_startDetectionTime+1+_d))
		throw RuntimeException("SMCAlgorithm::Run: Not enough observations.");
	
	// memory is reserved
	for(int iParticle=0;iParticle<_nParticles;iParticle++)
	{
		_estimatedChannelMatrices[iParticle] = new tMatrix[_endDetectionTime];
		_detectedSymbols[iParticle] = *(new tMatrix(_N,_endDetectionTime));
	}

// 	_estimatedChannelMatrices[0][0] = *(new tMatrix(2,2));
// 	_estimatedChannelMatrices[0][0](0,0) = 1;_estimatedChannelMatrices[0][0](0,1) = 2;	
// 	cout << "No llamo a process" << endl;

	this->Process(observations,noiseVariances);
}

void SMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
	if(observations.rows()!=_L || trainingSequence.rows()!=_N)
		throw RuntimeException("Run: Observations matrix or training sequence dimensions are wrong.");

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

void SMCAlgorithm::Resampling(int endResamplingTime)
{
	if(_resamplingCriterion.ResamplingNeeded(_weights))
	{
		vector<int> indexes = StatUtil::Discrete_rnd(_nParticles,_weights);
		_resamplingAlgorithm.Resampling(&_estimatedChannelMatrices,&_detectedSymbols,&_particlesChannelMatrixEstimators,indexes,_nParticles,_startDetectionTime,endResamplingTime,_endDetectionTime);
		_weights = 1.0/(double)_nParticles;
// 		cout << "Resampling..." << endl;
	}
}

double SMCAlgorithm::SER(tMatrix symbols)
{
	int windowSize = symbols.cols();
	int nDetectedVectors = _detectedSymbols[0].cols();

	if(windowSize>nDetectedVectors)
		throw RuntimeException("SMCAlgorithm::SER: more symbol vectors passed than detected.");
	
	// best particle is chosen
	int iBestParticle;
	Util::Max(_weights,iBestParticle);

	int nErrors = 0;
	int windowStart = nDetectedVectors - windowSize;
	int j;

	for(int i=windowStart;i<nDetectedVectors;i++)
	{
		j=0;
		while(j<_N)
		{
			if(_detectedSymbols[iBestParticle](j,i)!=symbols(j,i-windowStart))
				nErrors++;
			j++;
		}
	}
	return ((double)nErrors)/(double)(windowSize*_N);
}

