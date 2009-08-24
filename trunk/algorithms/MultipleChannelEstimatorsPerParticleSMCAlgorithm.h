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
    int _d;
    int _startDetectionTime;
    tRange _allSymbolsRows;

    double _channelUniqueMean, _channelUniqueVariance;
    vector<tMatrix> _channelMatrixMeans;
    vector<tMatrix> _channelMatrixVariances;

    vector<tMatrix> _channelMeanVectors;
    vector<tMatrix> _channelCovariances;
    
    bool _randomParticlesInitilization;

    virtual ParticleFilter* getParticleFilterPointer() = 0;
    virtual void initializeParticles() = 0;
    virtual void process(const tMatrix &observations,vector<double> noiseVariances) = 0;
    virtual int iBestChannelOrder(int iBestParticle) = 0;

    virtual void beforeInitializingParticles(const tMatrix &observations,vector<double> &noiseVariances,const tMatrix &trainingSequence) {}
    virtual void updateParticleChannelOrderEstimators(Particle *particle,const tMatrix &observations,const std::vector<std::vector<tMatrix> > &channelMatrices,vector<double> &noiseVariances,const tMatrix &sequenceToProcess) {}
public:
    MultipleChannelEstimatorsPerParticleSMCAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation,int smoothingLag,int nParticles,ResamplingAlgorithm *resamplingAlgorithm);

    void run(tMatrix observations,vector<double> noiseVariances);
    void run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);
    tMatrix getDetectedSymbolVectors();
    vector<tMatrix> getEstimatedChannelMatrices();

};

#endif
