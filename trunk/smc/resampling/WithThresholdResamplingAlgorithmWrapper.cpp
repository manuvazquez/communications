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
	#ifdef DEBUG
		cout << "Constructor" << endl;
	#endif
}

WithThresholdResamplingAlgorithmWrapper* WithThresholdResamplingAlgorithmWrapper::Clone() const
{
	return new WithThresholdResamplingAlgorithmWrapper(*this);
}

WithThresholdResamplingAlgorithmWrapper::~WithThresholdResamplingAlgorithmWrapper()
{
	delete _realResamplingAlgorithm;
}

WithThresholdResamplingAlgorithmWrapper::WithThresholdResamplingAlgorithmWrapper(const WithThresholdResamplingAlgorithmWrapper& withThresholdResamplingAlgorithmWrapper):ResamplingAlgorithm(withThresholdResamplingAlgorithmWrapper),_threshold(withThresholdResamplingAlgorithmWrapper._threshold),_realResamplingAlgorithm(withThresholdResamplingAlgorithmWrapper._realResamplingAlgorithm->Clone())
{
}

tVector WithThresholdResamplingAlgorithmWrapper::FlattenWeights(const tVector &weights, double threshold) const
{
	double remainingWeight = 0.0;
	int nParticlesOverThreshold = 0;

	tVector newWeights = weights;

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

// void WithThresholdResamplingAlgorithmWrapper::Resample(ParticleFilter* particleFilter, const tVector& weights)
// {
// 	double remainingWeight = 0.0;
// 	int nParticlesOverThreshold = 0;
//
// 	tVector newWeights = weights;
//
// 	for(int iParticle=0;iParticle<particleFilter->Nparticles();iParticle++)
// 		if(newWeights(iParticle)>=_threshold)
// 		{
// 			nParticlesOverThreshold++;
// 			remainingWeight += (newWeights(iParticle) - _threshold);
// 			newWeights(iParticle) = _threshold;
// 		}
//
// 	if(nParticlesOverThreshold > 0)
// 	{
// 		double weightToAddToEachParticle = remainingWeight/double(particleFilter->Nparticles() - nParticlesOverThreshold);
//
// 		for(int iParticle=0;iParticle<particleFilter->Nparticles();iParticle++)
// 			// if it wasn't a "<=" the obtained weights wouldn't be normalized
// 			if(newWeights(iParticle) < _threshold)
// 				newWeights(iParticle) += weightToAddToEachParticle;
// 	}
//
// 	#ifdef DEBUG
// 		cout << "Los pesos originales" << endl << weights;
// 		cout << "Los pesos modificados" << endl << newWeights;
// 		cout << "Su suma: " << Util::Sum(newWeights) << endl;
// 	#endif
//
// 	_realResamplingAlgorithm->Resample(particleFilter,newWeights);
//
// }

