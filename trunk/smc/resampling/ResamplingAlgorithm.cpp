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

#include "ResamplingAlgorithm.h"

// #define PRINT_RESAMPLED_INFO

bool ResamplingAlgorithm::ResampleWhenNecessary(ParticleFilter *particleFilter)
{
    tVector weigths = particleFilter->GetWeightsVector();

    if(_resamplingCriterion.resamplingNeeded(weigths))
    {
        vector<int> indexes = ObtainIndexes(particleFilter->Capacity(),weigths);
        particleFilter->keepParticles(indexes);

#ifdef PRINT_RESAMPLED_INFO
        cout << "particles that are going to be resampled..." << endl;
        vector<int> occurrences;
        vector<int> times;
        Util::howManyTimes(indexes,occurrences,times);
        cout << "occurrences (" << occurrences.size() << ")\t / times" << endl;
        for(uint i=0;i<occurrences.size();i++)
            cout << occurrences[i] << "\t" << times[i] << endl;
#endif
        return true;
    }
    return false;
}
