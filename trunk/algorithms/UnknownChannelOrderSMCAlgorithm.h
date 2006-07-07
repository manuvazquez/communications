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
#ifndef UNKNOWNCHANNELORDERSMCALGORITHM_H
#define UNKNOWNCHANNELORDERSMCALGORITHM_H

#include <UnknownChannelOrderAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <ParticleFilter.h>
#include <StdResamplingAlgorithm.h>
#include <ResamplingCriterion.h>
#include <ParticleWithChannelEstimationAndChannelOrder.h>

class UnknownChannelOrderSMCAlgorithm : public UnknownChannelOrderAlgorithm
{
protected:
	ParticleFilter _particleFilter;
    int _d,_startDetectionObservation,_startDetectionSymbolVector;
	double *_channelOrderWeightsSum;
    tRange _allSymbolsRows;
	vector<int> _nParticlesPerChannelOrder;

	virtual void InitializeParticles();
    virtual void Process(const tMatrix &observations,vector<double> noiseVariances) = 0;
	virtual void Resampling();
	void NormalizeParticleGroups();
public:
    UnknownChannelOrderSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int firstObservationIndex,int smoothingLag,int nParticles,ResamplingCriterion resamplingCriterion,StdResamplingAlgorithm resamplingAlgorithm);

    ~UnknownChannelOrderSMCAlgorithm();

	void Run(tMatrix observations,vector<double> noiseVariances);
    void Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);
	tMatrix GetDetectedSymbolVectors() {return tMatrix(0,0);}
	vector<tMatrix> GetEstimatedChannelMatrices() {return vector<tMatrix>(2);}
};

#endif
