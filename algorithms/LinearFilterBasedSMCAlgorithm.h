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
#include <ARchannel.h>
#include <ChannelDependentNoise.h>
#include <ParticleWithChannelEstimationAndLinearDetection.h>

class LinearFilterBasedSMCAlgorithm : public SMCAlgorithm
{
public:
    LinearFilterBasedSMCAlgorithm(std::string name, Alphabet alphabet,uint L,uint Nr,uint N, uint iLastSymbolVectorToBeDetected,uint m, ChannelMatrixEstimator *channelEstimator,LinearDetector *linearDetector,MatrixXd preamble, uint SMCsmoothingLag, uint nParticles, ResamplingAlgorithm *resamplingAlgorithm,const MatrixXd &channelMatrixMean, const MatrixXd &channelMatrixVariances,double ARcoefficient,double samplingVariance, double ARprocessVariance, bool substractContributionFromKnownSymbols=false);

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
    LinearFilterBasedSMCAlgorithm(std::string name, Alphabet alphabet,uint L,uint Nr,uint N, uint iLastSymbolVectorToBeDetected,uint m,MatrixXd preamble, uint SMCsmoothingLag, ParticleFilter *particleFilter, ResamplingAlgorithm *resamplingAlgorithm,double ARcoefficient,double samplingVariance, double ARprocessVariance, bool substractContributionFromKnownSymbols=false);

    ~LinearFilterBasedSMCAlgorithm();

protected:
    LinearDetector *_linearDetector;
    double _ARcoefficient,_samplingVariance,_ARprocessVariance;

    void initializeParticles();
    void process(const MatrixXd &observations, vector<double> noiseVariances);

    bool _substractContributionFromKnownSymbols;

    virtual void fillFirstEstimatedChannelMatrix(uint iParticle,MatrixXd &firstEstimatedChannelMatrix) const
    {
        // firstEstimatedChannelMatrix = _ARcoefficient * <lastEstimatedChannelMatrix> + randn(_nOutputs,_nInputsXchannelOrder)*_samplingVariance
        firstEstimatedChannelMatrix = _ARcoefficient*dynamic_cast<WithChannelEstimationParticleAddon *>(_particleFilter->getParticle(iParticle))->getChannelMatrixEstimator(_estimatorIndex)->lastEstimatedChannelMatrix() + StatUtil::randnMatrix(_nOutputs,_nInputsXchannelOrder,0.0,_samplingVariance);
        
    }

    virtual void beforeInitializingParticles(const MatrixXd &observations, const MatrixXd &trainingSequence);
};

#endif
