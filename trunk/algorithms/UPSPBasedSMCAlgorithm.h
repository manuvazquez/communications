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
#ifndef UPSPBASEDSMCALGORITHM_H
#define UPSPBASEDSMCALGORITHM_H

#include <MultipleChannelEstimatorsPerParticleSMCAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/
class UPSPBasedSMCAlgorithm : public MultipleChannelEstimatorsPerParticleSMCAlgorithm
{
protected:
	ParticleFilter *_particleFilter;
	double _ARcoefficient,_samplingVariance,_ARprocessVariance;

    int BestChannelOrderIndex(int iBestParticle);
public:
    UPSPBasedSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance);

    ~UPSPBasedSMCAlgorithm();

    virtual ParticleFilter* GetParticleFilterPointer() {return _particleFilter;}
    virtual void InitializeParticles();
    virtual void Process(const tMatrix& observations, vector< double > noiseVariances);

};

#endif