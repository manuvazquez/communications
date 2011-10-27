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
#ifndef ISIS_H
#define ISIS_H

#include <MultipleChannelEstimatorsPerParticleSMCAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <ParticleWithChannelEstimationAndChannelOrderAPP.h>
#include <KalmanEstimator.h>

class ISIS : public MultipleChannelEstimatorsPerParticleSMCAlgorithm
{
protected:
	ParticleFilter _particleFilter;

    virtual ParticleFilter* getParticleFilterPointer() {return &_particleFilter;}
    virtual void initializeParticles();
    
    virtual void process(const MatrixXd& observations, vector<double> noiseVariances);

	// it's never gonna be called
    uint iBestChannelOrder(uint iBestParticle) { return 0u;}
public:
    ISIS(string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, MatrixXd preamble, uint iFirstObservation, uint smoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm);

    vector<MatrixXd> getEstimatedChannelMatrices();
};

#endif
