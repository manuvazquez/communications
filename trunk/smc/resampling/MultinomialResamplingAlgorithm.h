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
#ifndef MULTINOMIALRESAMPLINGALGORITHM_H
#define MULTINOMIALRESAMPLINGALGORITHM_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <types.h>
#include <ResamplingAlgorithm.h>
// #include <ChannelMatrixEstimator.h>
// #include <LinearDetector.h>
// #include <ParticleFilter.h>
// #include <ParticleWithChannelEstimation.h>
// #include <exceptions.h>

class MultinomialResamplingAlgorithm : public ResamplingAlgorithm{
public:
    MultinomialResamplingAlgorithm(ResamplingCriterion resamplingCriterion);

	virtual MultinomialResamplingAlgorithm* Clone() const;

	void Resample(ParticleFilter *particleFilter,const tVector &weights);
// 	int ResampleWhenNecessary(ParticleFilter *particleFilter);
};

#endif
