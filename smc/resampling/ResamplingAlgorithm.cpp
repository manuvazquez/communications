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

bool ResamplingAlgorithm::resampleWhenNecessary(ParticleFilter *particleFilter)
{
    VectorXd weigths = particleFilter->getWeightsVector();

    if(_resamplingCriterion.resamplingNeeded(weigths))
    {
        vector<uint> indexes = obtainIndexes(particleFilter->capacity(),weigths);
        particleFilter->keepParticles(indexes);

#ifdef PRINT_RESAMPLED_INFO
        cout << "particles that are going to be resampled..." << endl;
        vector<uint> occurrences;
        vector<uint> times;
        Util::howManyTimes(indexes,occurrences,times);
        cout << "occurrences (" << occurrences.size() << ")\t / times" << endl;
        for(uint i=0;i<occurrences.size();i++)
            cout << occurrences[i] << "\t" << times[i] << endl;
        cout << "stroke to continue..." << endl;
        getchar();
#endif
        return true;
    }
    return false;
}

std::vector<uint> ResamplingAlgorithm::obtainIndexes(uint n,const VectorXd &weights, const std::vector<bool> &mask)
{
	if(static_cast<uint>(weights.size())!=mask.size())
	  throw RuntimeException("ResamplingAlgorithm::obtainIndexes: mask has not the same dimensions as the weights vector");

	uint nActiveElements = 0;
	for(uint i=0;i<mask.size();i++)
		nActiveElements += mask[i];
	
	// if no element is to be taken into account...
	if(nActiveElements==0)
		//...we directly return an empty vector
		return std::vector<uint>(0u);
	
	VectorXd filteredWeights(nActiveElements);
	std::vector<uint> mapping;
	
	double normConst = 0.0;
  
	nActiveElements = 0;
	for(uint i=0;i<weights.size();i++)
		if(mask[i])
		{
			filteredWeights(nActiveElements++) = weights(i);
			mapping.push_back(i);

			normConst += weights(i);
		}
	  
	filteredWeights /= normConst;

	std::vector<uint> relativeIndexes = obtainIndexes(n,filteredWeights);
	
// 	cout << "relative indexes" << endl;
// 	Util::print(relativeIndexes);
// 	cout << endl;
// 
// 	cout << "mapping" << endl;
// 	Util::print(mapping);
// 	cout << endl;

	std::vector<uint> res(relativeIndexes.size());
	for(uint i=0;i<relativeIndexes.size();i++)
		res[i] = mapping[relativeIndexes[i]];
	
	return res;
}
