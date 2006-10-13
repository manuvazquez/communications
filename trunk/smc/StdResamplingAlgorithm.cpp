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
#include "StdResamplingAlgorithm.h"

// #define DEBUG

using namespace std;

StdResamplingAlgorithm::StdResamplingAlgorithm(ResamplingCriterion resamplingCriterion):ResamplingAlgorithm(resamplingCriterion)
{
}

void StdResamplingAlgorithm::Resample(ParticleFilter *particleFilter)
{
    tVector weigths = particleFilter->GetWeightsVector();

    if(_resamplingCriterion.ResamplingNeeded(weigths))
    {
        vector<int> indexes = StatUtil::Discrete_rnd(particleFilter->Nparticles(),weigths);
        particleFilter->KeepParticles(indexes);

		#ifdef DEBUG
		for(int i=0;i<particleFilter->Nparticles();i++)
			cout << i << " <- " << indexes[i] << endl;
		#endif
    }
}
