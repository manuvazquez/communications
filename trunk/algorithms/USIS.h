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
#ifndef USIS_H
#define USIS_H

#include <MultipleChannelEstimatorsPerParticleSMCAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <string>
#include <sys/time.h>

#include <LinearDetector.h>
#include <ChannelOrderEstimator.h>
#include <ParticleFilter.h>
#include <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation.h>
#include <MIMOChannel.h>

class USIS : public MultipleChannelEstimatorsPerParticleSMCAlgorithm
{
public:
    USIS(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators,vector<LinearDetector *> linearDetectors, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, ChannelOrderEstimator * channelOrderEstimator, double ARcoefficient,double samplingVariance,double ARprocessVariance/*,const MIMOChannel &canal,const tMatrix &simbolos*/);

    ~USIS();

	#ifdef CHANNELORDERSAPP_SAVING
    	tMatrix GetChannelOrderAPPsAlongTime() { return _channelOrderAPPs(_rCandidateOrders,tRange(_preamble.cols(),_K-1));}
    #endif

    bool PerformsChannelOrderAPPEstimation() { return true;}

protected:
    vector<LinearDetector *> _linearDetectors;
    ChannelOrderEstimator *_channelOrderEstimator;
	ParticleFilter _particleFilter;
	double _ARcoefficient,_samplingVariance,_ARprocessVariance;
    tRange _rAllObservationRows;
    #ifdef CHANNELORDERSAPP_SAVING
		tMatrix _channelOrderAPPs;
	#endif
	tRange _rCandidateOrders;
	bool _processDoneExternally;

    virtual ParticleFilter* GetParticleFilterPointer() {return &_particleFilter;}
    vector<vector<tMatrix> > ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence);
    virtual void InitializeParticles();
    virtual void Process(const tMatrix& observations, vector< double > noiseVariances);

	virtual void BeforeResamplingProcess(int iProcessedObservation, const tMatrix& observations, const vector< double > &noiseVariances) {}

	vector<tMatrix> GetEstimatedChannelMatrices();
};

#endif
