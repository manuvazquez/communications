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

#include <ParticleFilter.h>
#include <ResamplingCriterion.h>
#include <ResamplingAlgorithm.h>
#include <StatUtil.h>
#include <ParticleWithChannelEstimation.h>
#include <math.h>

/**
It implements some parts of a typical (default) SMC algorithm. It assumes that the particles have a channel estimator (derive from the class "ParticleWithChannelEstimation) and thus, methods "initializeParticles" or "run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)" must be redefined in a subclass implementing an algorithm which requires another type of particle.

    @author Manu <manu@rustneversleeps>
*/
class SMCAlgorithm : public KnownChannelOrderAlgorithm
{
protected:
    ParticleFilter *_particleFilter;
    bool _particleFilterNeedToBeDeleted;
    ResamplingAlgorithm *_resamplingAlgorithm;
    uint _d,_startDetectionTime;

    // a particle contains a vector of channel estimators (and possibly linear detectors)
    uint _estimatorIndex; //! it indicates which of the all the estimator that each particle contain is interesting at every moment

    VectorXd _channelMean;
    MatrixXd _channelCovariance;

    virtual void initializeParticles();
    virtual void process(const MatrixXd &observations,vector<double> noiseVariances) = 0;
    
    /**
     *    Computes the smoothed likelihood, i.e., the product of the likelihoods of the observations involved in the smoothing from @param iObservationToBeProcessed to @param iObservationToBeProcessed + d (the smoothing lag)
     * @param channelMatrices
     * @param involvedSymbolVectors
     * @param particle
     * @param iObservationToBeProcessed
     * @param observations
     * @param noiseVariances
     * @return
     */
    double smoothedLikelihood(const vector<MatrixXd> &channelMatrices,const MatrixXd &involvedSymbolVectors,uint iObservationToBeProcessed,const MatrixXd &observations,const vector<double> &noiseVariances);

    bool _randomParticlesInitilization;

    virtual void beforeInitializingParticles(const MatrixXd &observations, const MatrixXd &trainingSequence) {}

public:
    SMCAlgorithm(string name, Alphabet alphabet,uint L,uint Nr,uint N, uint iLastSymbolVectorToBeDetected,uint m, ChannelMatrixEstimator *channelEstimator, MatrixXd preamble,uint smoothingLag,uint nParticles,ResamplingAlgorithm *resamplingAlgorithm, const MatrixXd &channelMatrixMean, const MatrixXd &channelMatrixVariances);

    /**
     * Constructor for allowing the algorithm to operate over a already constructed particle filter
     * @param name
     * @param alphabet
     * @param L
     * @param N
     * @param iLastSymbolVectorToBeDetected
     * @param m
     * @param channelEstimator
     * @param preamble
     * @param smoothingLag
     * @param particleFilter
     * @param resamplingAlgorithm
     */
    SMCAlgorithm(string name, Alphabet alphabet,uint L,uint Nr,uint N, uint iLastSymbolVectorToBeDetected,uint m, MatrixXd preamble,uint smoothingLag,ParticleFilter *particleFilter,ResamplingAlgorithm *resamplingAlgorithm);

    ~SMCAlgorithm();

    void setEstimatorIndex(uint n);
    
    virtual void run(MatrixXd observations,vector<double> noiseVariances);
    virtual void runFrom(uint n,MatrixXd observations,vector<double> noiseVariances);    
    virtual void run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence);

    virtual MatrixXd getDetectedSymbolVectors();
    
    virtual vector<MatrixXd> getEstimatedChannelMatrices();    
};

#endif
