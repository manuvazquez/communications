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
#ifndef UCO_SIS_H
#define UCO_SIS_H

#include <MultipleChannelEstimatorsPerParticleSMCAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <LinearDetector.h>
#include <ParticleFilter.h>
#include <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP.h>
#include <MIMOChannel.h>

class UCO_SIS : public MultipleChannelEstimatorsPerParticleSMCAlgorithm
{
public:
    UCO_SIS(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators,vector<LinearDetector *> linearDetectors, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance/*,const MIMOChannel &canal,const tMatrix &simbolos*/);

    ~UCO_SIS();

protected:
    vector<LinearDetector *> _linearDetectors;
	ParticleFilter _particleFilter;
	double _ARcoefficient,_samplingVariance,_ARprocessVariance;
    tRange _rAllObservationRows;

// 	const MIMOChannel &_canal;
// 	const tMatrix &_simbolos;

    virtual ParticleFilter* GetParticleFilterPointer() {return &_particleFilter;}
    vector<vector<tMatrix> > ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence);
    virtual void InitializeParticles();
    virtual void Process(const tMatrix& observations, vector< double > noiseVariances);
};

#endif
