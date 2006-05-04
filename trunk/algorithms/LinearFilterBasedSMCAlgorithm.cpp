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

LinearFilterBasedSMCAlgorithm::LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet, ChannelMatrixEstimator& channelEstimator,LinearDetector &linearDetector,tMatrix preamble, int smoothingLag, int nParticles, ResamplingCriterion resamplingCriterion, StdResamplingAlgorithm resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance): SMCAlgorithm(name, alphabet, channelEstimator, preamble, smoothingLag, nParticles, resamplingCriterion, resamplingAlgorithm),_particlesLinearDetectors(new LinearDetector*[_nParticles]),_linearDetector(&linearDetector),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance)
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

void LinearFilterBasedSMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
	for(int iParticle=0;iParticle<_nParticles;iParticle++)
	{
		_particlesLinearDetectors[iParticle] = _linearDetector->Clone();
	}

	SMCAlgorithm::Run(observations,noiseVariances);
}

void LinearFilterBasedSMCAlgorithm::Process(tMatrix observations, vector< double > noiseVariances)
{
	cout << "Processing in LinearFilterBasedSMCAlgorithm..." << endl;
	int iParticle,iSmoothing;
	vector<tMatrix> matricesToStack(_d+1,tMatrix(_L,_Nm));
	tRange allObservationsRows(0,_L-1);

	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_endDetectionTime;iObservationToBeProcessed++)
	{
		// observation matrix columns that are involved in the smoothing
		tRange smoothingRange(iObservationToBeProcessed,iObservationToBeProcessed+_d);

		// the stacked observations vector
		tVector stackedObservations = Util::ToVector(observations(allObservationsRows,smoothingRange),columnwise);		

		for(iParticle=0;iParticle<_nParticles;iParticle++)
		{
// 			// predicted channel matrices are stored in a vector in order to stack them

			// matricesToStack[0] = _ARcoefficient * <lastEstimatedChannelMatrix> + randn(_L,_Nm)*_samplingVariance
			Util::Add((_particlesChannelMatrixEstimators[iParticle])->LastEstimatedChannelMatrix(),StatUtil::RandnMatrix(_L,_Nm,0.0,_samplingVariance),matricesToStack[0],_ARcoefficient,1.0);

			for(iSmoothing=1;iSmoothing<=_d;iSmoothing++)
			{
				// matricesToStack[iSmoothing] = _ARcoefficient * matricesToStack[iSmoothing-1] + rand(_L,_Nm)*_ARprocessVariance
				Util::Add(matricesToStack[iSmoothing-1],StatUtil::RandnMatrix(_L,_Nm,0.0,_ARprocessVariance),matricesToStack[iSmoothing],_ARcoefficient,1.0);
			}

			// matrices are stacked to give
			tMatrix stackedChannelMatrix = HsToStackedH(matricesToStack);

			// the estimated stacked channel matrix is used to obtain soft estimations
			// of the transmitted symbols
			tVector softEstimations =  (_particlesLinearDetectors[iParticle])->Detect(stackedObservations,stackedChannelMatrix);
		}
	}
}

vector<tMatrix> LinearFilterBasedSMCAlgorithm::ProcessTrainingSequence(tMatrix observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
	int lengthSequenceToProcess = _preamble.cols() + trainingSequence.cols();
	tRange allObservationRows(0,_L-1);

	for(int i=_m-1;i<lengthSequenceToProcess;i++)
	{
		tRange smoothingRange(i,i+_d);
		tVector stackedObservationsVector = Util::ToVector(observations(allObservationRows,smoothingRange),columnwise);
		_linearDetector->StateStep(stackedObservationsVector);
	}

	// the resultant linear detector is copied into each particle
	for(int iParticle=0;iParticle<_nParticles;iParticle++)
		_particlesLinearDetectors[iParticle] = _linearDetector->Clone();

	return SMCAlgorithm::ProcessTrainingSequence(observations,noiseVariances,trainingSequence);
}
