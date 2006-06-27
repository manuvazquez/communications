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
#ifndef LINEARFILTERBASEDSMCALGORITHM_H
#define LINEARFILTERBASEDSMCALGORITHM_H

#include <SMCAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <LinearDetector.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <ARchannel.h>
#include <ChannelDependentNoise.h>
#include <ParticleWithChannelEstimationAndLinearDetection.h>

class LinearFilterBasedSMCAlgorithm : public SMCAlgorithm
{
public:
    LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int N, int K, ChannelMatrixEstimator *channelEstimator,LinearDetector *linearDetector,tMatrix preamble, int smoothingLag, int nParticles, ResamplingCriterion resamplingCriterion, StdResamplingAlgorithm resamplingAlgorithm,double ARcoefficient,double samplingVariance, double ARprocessVariance);

    ~LinearFilterBasedSMCAlgorithm();

// 	using SMCAlgorithm::Run;
// 	void Run(tMatrix observations,vector<double> noiseVariances);

protected:
// 	LinearDetector **_particlesLinearDetectors;
	LinearDetector *_linearDetector;
	double _ARcoefficient,_samplingVariance,_ARprocessVariance;

	void InitializeParticles();
    void Process(const tMatrix &observations, vector< double > noiseVariances);
	vector<tMatrix> ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence);
// 	void Resampling(int endResamplingTime);
};

#endif
