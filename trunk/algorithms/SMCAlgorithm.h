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
It implements some parts of a typical (default) SMC algorithm. It assumes that the particles have a channel estimator (derive from the class "ParticleWithChannelEstimation) and thus, methods "InitializeParticles" or "Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)" must be redefined in a subclass implementing an algorithm which requires another type of particle.

    @author Manu <manu@rustneversleeps>
*/
class SMCAlgorithm : public KnownChannelOrderAlgorithm
{
protected:
    ParticleFilter *_particleFilter;
    bool _particleFilterNeedToBeDeleted;
    ResamplingAlgorithm *_resamplingAlgorithm;
    int _d,_startDetectionTime;
    tRange _allSymbolsRows;

    // a particle contains a vector of channel estimators (and possibly linear detectors)
    int _estimatorIndex; //! it indicates which of the all the estimator that each particle contain is interesting at every moment

    tMatrix _channelMatrixMean,_channelMatrixVariances;

    virtual void InitializeParticles();
    virtual void Process(const tMatrix &observations,vector<double> noiseVariances) = 0;
    
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
    double smoothedLikelihood(const vector<tMatrix> &channelMatrices,const tMatrix &involvedSymbolVectors,ParticleWithChannelEstimation *particle,int iObservationToBeProcessed,const tMatrix &observations,const vector<double> &noiseVariances);

    const MIMOChannel *_channel;
    const tMatrix *_symbols;

    bool _randomParticlesInitilization;

    virtual void BeforeInitializingParticles(const tMatrix &observations, const tMatrix &trainingSequence) {}

public:
    SMCAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m, ChannelMatrixEstimator *channelEstimator, tMatrix preamble,int smoothingLag,int nParticles,ResamplingAlgorithm *resamplingAlgorithm, const tMatrix &channelMatrixMean, const tMatrix &channelMatrixVariances);

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
    SMCAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m, tMatrix preamble,int smoothingLag,ParticleFilter *particleFilter,ResamplingAlgorithm *resamplingAlgorithm);

    ~SMCAlgorithm();

    void SetEstimatorIndex(int n);

    void Run(tMatrix observations,vector<double> noiseVariances);
    void Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);

    void RunFrom(int n,tMatrix observations,vector<double> noiseVariances);

    tMatrix getDetectedSymbolVectors();
    vector<tMatrix> GetEstimatedChannelMatrices();
};

#endif
