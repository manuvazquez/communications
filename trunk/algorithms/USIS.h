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
    USIS(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators,vector<LinearDetector *> linearDetectors, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, ChannelOrderEstimator * channelOrderEstimator, double ARcoefficient,double samplingVariance,double ARprocessVariance);

    ~USIS();

protected:
    vector<LinearDetector *> _linearDetectors;
    ChannelOrderEstimator *_channelOrderEstimator;
    ParticleFilter _particleFilter;
    double _ARcoefficient,_samplingVariance,_ARprocessVariance;
    bool _processDoneExternally;

    virtual ParticleFilter* getParticleFilterPointer() {return &_particleFilter;}
    virtual void initializeParticles();
    
    virtual void process(const tMatrix& observations, vector< double > noiseVariances)
    {
        process(Util::lapack2eigen(observations),noiseVariances);
    }
    virtual void process(const MatrixXd& observations, vector< double > noiseVariances);

    virtual void beforeResamplingProcess(int iProcessedObservation, const tMatrix& observations, const vector<double> &noiseVariances) {}
    virtual void beforeResamplingProcess(int iProcessedObservation, const MatrixXd& observations, const vector<double> &noiseVariances) {}

    int iBestChannelOrder(int iBestParticle);

    virtual void beforeInitializingParticles(const tMatrix &observations,vector<double> &noiseVariances,const tMatrix &trainingSequence)
    {
        beforeInitializingParticles(Util::lapack2eigen(observations),noiseVariances,Util::lapack2eigen(trainingSequence));
    }
    virtual void beforeInitializingParticles(const MatrixXd &observations, vector<double> &noiseVariances, const MatrixXd &trainingSequence);
    
    virtual void updateParticleChannelOrderEstimators(Particle *particle,const tMatrix &observations,const std::vector<std::vector<tMatrix> > &channelMatrices,vector<double> &noiseVariances,const tMatrix &sequenceToProcess)
    {
        updateParticleChannelOrderEstimators(particle,Util::lapack2eigen(observations),Util::lapack2eigen(channelMatrices),noiseVariances,Util::lapack2eigen(sequenceToProcess));
    }
    virtual void updateParticleChannelOrderEstimators(Particle *particle,const MatrixXd &observations,const std::vector<std::vector<MatrixXd> > &channelMatrices,vector<double> &noiseVariances,const MatrixXd &sequenceToProcess);
};

#endif
