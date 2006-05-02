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
#include "LinearFilterBasedSMCAlgorithm.h"

LinearFilterBasedSMCAlgorithm::LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet, ChannelMatrixEstimator& channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingCriterion resamplingCriterion, StdResamplingAlgorithm resamplingAlgorithm,LinearDetector &linearDetector): SMCAlgorithm(name, alphabet, channelEstimator, preamble, smoothingLag, nParticles, resamplingCriterion, resamplingAlgorithm),_particlesLinearDetectors(new LinearDetector*[_nParticles]),_linearDetector(&linearDetector)
{
	for(int i=0;i<_nParticles;i++)
	{
		_particlesLinearDetectors[i] = NULL;
	}
}


LinearFilterBasedSMCAlgorithm::~LinearFilterBasedSMCAlgorithm()
{
	for(int i=0;i<_nParticles;i++)
	{
		delete _particlesLinearDetectors[i];
	}
	delete[] _particlesLinearDetectors;
}


void LinearFilterBasedSMCAlgorithm::Process(tMatrix observations, vector< double > noiseVariances)
{
}

vector<tMatrix> LinearFilterBasedSMCAlgorithm::ProcessTrainingSequence(tMatrix observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
	int lengthSequenceToProcess = _preamble.cols() + trainingSequence.cols();
	tRange allObservationsRows(0,_L-1);

	// the channel matrix received by the method Detect doesn't affect the state of the detector
	tMatrix nullChannelMatrix(_L*(_d+1),_N*(_m+_d));

	for(int i=_m-1;i<lengthSequenceToProcess;i++)
	{
		tRange smoothingRange(i,i+_d);
// 		tVector stackedObservationsVector = Util::ToVector(observations(allObservationsRows,smoothingRange),columnwise);
// 		_linearDetector->Detect(stackedObservationsVector,nullChannelMatrix);
// 		_linearDetector->Detect(observations.col(i),
	}
}
