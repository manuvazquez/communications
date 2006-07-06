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
#ifndef RESAMPLINGCRITERION_H
#define RESAMPLINGCRITERION_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <types.h>
#include "smcExceptions.h"
#include <ParticleWithChannelEstimation.h>

class ResamplingCriterion{
private:
	double _resamplingRatio;
public:
    ResamplingCriterion(double resamplingRatio);

	bool ResamplingNeeded(tVector weights,std::vector<int> indexes);

	bool ResamplingNeeded(tVector weights)
	{
		int nParticles = weights.size();
		std::vector<int> indexes(nParticles);
		for(int i=0;i<nParticles;i++)
			indexes[i] = i;
		return ResamplingNeeded(weights,indexes);
	}

// 	bool ResamplingNeeded(ParticleWithChannelEstimation **particles,std::vector<int> indexes);
// 
// 	bool ResamplingNeeded(ParticleWithChannelEstimation **particles,int nParticles)
// 	{
// 		std::vector<int> indexes(nParticles);
// 		for(int i=0;i<nParticles;i++)
// 			indexes[i] = i;
// 		return ResamplingNeeded(particles,indexes);
// 	}

};

#endif
