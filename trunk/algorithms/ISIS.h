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
    int iBestChannelOrder(int iBestParticle) { return 0;}
public:
    ISIS(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm);

    vector<tMatrix> getEstimatedChannelMatrices()
    {
        return Util::eigen2lapack(getEstimatedChannelMatrices_eigen());
    }
    vector<MatrixXd> getEstimatedChannelMatrices_eigen();
};

#endif
