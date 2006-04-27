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
#ifndef SMCALGORITHM_H
#define SMCALGORITHM_H

#include <KnownChannelOrderAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <ResamplingCriterion.h>
#include <StdResamplingAlgorithm.h>
#include <StatUtil.h>

class SMCAlgorithm : public KnownChannelOrderAlgorithm
{
protected:
	int _d, _nParticles, _startDetectionTime, _endDetectionTime;
	ResamplingCriterion _resamplingCriterion;
	StdResamplingAlgorithm _resamplingAlgorithm;
// 	vector<tMatrix>  *estimatedChannelMatrices;
	tMatrix  **_estimatedChannelMatrices;
	tMatrix **_detectedSymbols;
	ChannelMatrixEstimator **_particlesChannelMatrixEstimators;
	tVector _weights;
	tRange _allSymbolsRows;
	bool _reservedMemory;

	virtual void Resampling(int endResamplingTime);
	virtual void Process(tMatrix observations,vector<double> noiseVariances) = 0;
public:
    SMCAlgorithm(string name, Alphabet alphabet, ChannelMatrixEstimator& channelEstimator, tMatrix preamble,int smoothingLag,int nParticles,ResamplingCriterion resamplingCriterion,StdResamplingAlgorithm resamplingAlgorithm);
	~SMCAlgorithm();
	
	void Run(tMatrix observations,vector<double> noiseVariances);
	void Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);
	double SER(tMatrix symbols);
};

#endif