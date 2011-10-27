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
#ifndef LINEARFILTERBASEDMKFALGORITHM_H
#define LINEARFILTERBASEDMKFALGORITHM_H

#include <LinearFilterBasedSMCAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <KalmanEstimator.h>

class LinearFilterBasedMKFAlgorithm : public LinearFilterBasedSMCAlgorithm
{
public:
    LinearFilterBasedMKFAlgorithm(string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, KalmanEstimator* channelEstimator, LinearDetector* linearDetector, MatrixXd preamble, uint backwardsSmoothingLag, uint smoothingLag, int forwardSmoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm, const MatrixXd& channelMatrixMean, const MatrixXd& channelMatrixVariances, double ARcoefficient, double samplingVariance, double ARprocessVariance, bool substractContributionFromKnownSymbols=false);

protected:
    virtual void fillFirstEstimatedChannelMatrix(int iParticle, MatrixXd& firstEstimatedChannelMatrix) const
    {
        firstEstimatedChannelMatrix = (dynamic_cast<KalmanEstimator *> (dynamic_cast<WithChannelEstimationParticleAddon *>(_particleFilter->getParticle(iParticle))->getChannelMatrixEstimator(_estimatorIndex)))->sampleFromPredictive();
    }
};

#endif
