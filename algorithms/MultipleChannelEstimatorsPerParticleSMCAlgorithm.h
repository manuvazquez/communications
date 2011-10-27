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
#ifndef MULTIPLECHANNELESTIMATORSPERPARTICLESMCALGORITHM_H
#define MULTIPLECHANNELESTIMATORSPERPARTICLESMCALGORITHM_H

#include <UnknownChannelOrderAlgorithm.h>

/**
    @author Manu <manu@rustneversleeps>
*/

#include <ResamplingAlgorithm.h>

class MultipleChannelEstimatorsPerParticleSMCAlgorithm : public UnknownChannelOrderAlgorithm
{
protected:
    ResamplingAlgorithm *_resamplingAlgorithm;
    uint _d;
    uint _startDetectionTime;

	// mean and variance which will server to initialize the first channel matrix estimations
    double _channelUniqueMean, _channelUniqueVariance;

	// a vector of matrices, each one representing the mean of a channel estimation for one channel order
    vector<MatrixXd> _channelMeanVectors;
	
	// a vector of matrices, each one representing the covariance of a channel estimation for one channel order
    vector<MatrixXd> _channelCovariances;
    
    bool _randomParticlesInitilization;

    virtual ParticleFilter* getParticleFilterPointer() = 0;
    virtual void initializeParticles() = 0;
    virtual void process(const MatrixXd &observations,vector<double> noiseVariances) = 0; // eigen
    virtual int iBestChannelOrder(int iBestParticle) = 0;

    virtual void beforeInitializingParticles(const MatrixXd &observations,vector<double> &noiseVariances,const MatrixXd &trainingSequence) {}
    virtual void updateParticleChannelOrderEstimators(Particle *particle,const MatrixXd &observations,const std::vector<std::vector<MatrixXd> > &channelMatrices,vector<double> &noiseVariances,const MatrixXd &sequenceToProcess) {}
public:
    MultipleChannelEstimatorsPerParticleSMCAlgorithm(string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, MatrixXd preamble, uint iFirstObservation,uint smoothingLag,uint nParticles,ResamplingAlgorithm *resamplingAlgorithm);

    virtual void run(MatrixXd observations,vector<double> noiseVariances);
    virtual void run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence);

    virtual MatrixXd getDetectedSymbolVectors();
    
    virtual vector<MatrixXd> getEstimatedChannelMatrices();  

};

#endif
