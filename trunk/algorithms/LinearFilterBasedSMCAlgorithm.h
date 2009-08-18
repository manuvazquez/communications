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
    LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m, ChannelMatrixEstimator *channelEstimator,LinearDetector *linearDetector,tMatrix preamble, int backwardsSmoothingLag, int SMCsmoothingLag, int forwardSmoothingLag, int nParticles, ResamplingAlgorithm *resamplingAlgorithm,const tMatrix &channelMatrixMean, const tMatrix &channelMatrixVariances,double ARcoefficient,double samplingVariance, double ARprocessVariance, bool substractContributionFromKnownSymbols=false);

    /**
     * Constructor for allowing the algorithm to operate over an already constructed particle filter
     * @param name
     * @param alphabet
     * @param L
     * @param N
     * @param iLastSymbolVectorToBeDetected
     * @param m
     * @param preamble
     * @param SMCsmoothingLag
     * @param nParticles
     * @param resamplingAlgorithm
     * @param ARcoefficient
     * @param samplingVariance
     * @param ARprocessVariance
     */
    LinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m,tMatrix preamble, int SMCsmoothingLag, ParticleFilter *particleFilter, ResamplingAlgorithm *resamplingAlgorithm,double ARcoefficient,double samplingVariance, double ARprocessVariance, bool substractContributionFromKnownSymbols=false);

    ~LinearFilterBasedSMCAlgorithm();

protected:
    LinearDetector *_linearDetector;
    double _ARcoefficient,_samplingVariance,_ARprocessVariance;
    int _c,_e;

    void initializeParticles();
    void process(const tMatrix &observations, vector< double > noiseVariances);

    bool _substractContributionFromKnownSymbols;

    virtual void fillFirstEstimatedChannelMatrix(int iParticle,tMatrix &firstEstimatedChannelMatrix) const
    {
        // firstEstimatedChannelMatrix = _ARcoefficient * <lastEstimatedChannelMatrix> + randn(_nOutputs,_nInputsXchannelOrder)*_samplingVariance
        Util::add(dynamic_cast<WithChannelEstimationParticleAddon *>(_particleFilter->getParticle(iParticle))->getChannelMatrixEstimator(_estimatorIndex)->lastEstimatedChannelMatrix(),StatUtil::randnMatrix(_nOutputs,_nInputsXchannelOrder,0.0,_samplingVariance),firstEstimatedChannelMatrix,_ARcoefficient,1.0);
    }

    virtual void beforeInitializingParticles(const tMatrix &observations, const tMatrix &trainingSequence);
};

#endif
