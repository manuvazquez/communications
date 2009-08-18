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
#include "WithThresholdResamplingAlgorithmWrapper.h"

// #define DEBUG

WithThresholdResamplingAlgorithmWrapper::WithThresholdResamplingAlgorithmWrapper(ResamplingAlgorithm *resamplingAlgorithm,double threshold): ResamplingAlgorithm(resamplingAlgorithm->GetResamplingCriterion()),_threshold(threshold),_realResamplingAlgorithm(resamplingAlgorithm)
{
}

WithThresholdResamplingAlgorithmWrapper* WithThresholdResamplingAlgorithmWrapper::clone() const
{
	return new WithThresholdResamplingAlgorithmWrapper(*this);
}

WithThresholdResamplingAlgorithmWrapper::~WithThresholdResamplingAlgorithmWrapper()
{
	delete _realResamplingAlgorithm;
}

WithThresholdResamplingAlgorithmWrapper::WithThresholdResamplingAlgorithmWrapper(const WithThresholdResamplingAlgorithmWrapper& withThresholdResamplingAlgorithmWrapper):ResamplingAlgorithm(withThresholdResamplingAlgorithmWrapper),_threshold(withThresholdResamplingAlgorithmWrapper._threshold),_realResamplingAlgorithm(withThresholdResamplingAlgorithmWrapper._realResamplingAlgorithm->clone())
{
}

// eigen
VectorXd WithThresholdResamplingAlgorithmWrapper::FlattenWeights(const VectorXd &weights, double threshold) const
{
    double remainingWeight = 0.0;
    int nParticlesOverThreshold = 0;

    VectorXd newWeights = weights;

    for(int iWeight=0;iWeight<weights.size();iWeight++)
        if(newWeights(iWeight)>=threshold)
        {
            nParticlesOverThreshold++;
            remainingWeight += (newWeights(iWeight) - threshold);
            newWeights(iWeight) = threshold;
        }

    if(nParticlesOverThreshold > 0)
    {
        double weightToAddToEachParticle = remainingWeight/double(weights.size() - nParticlesOverThreshold);

        for(int iWeight=0;iWeight<weights.size();iWeight++)
            if(newWeights(iWeight) < threshold)
                newWeights(iWeight) += weightToAddToEachParticle;
    }

    return newWeights;
}
