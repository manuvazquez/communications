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
#include "ResidualResamplingAlgorithm.h"

// #define DEBUG

ResidualResamplingAlgorithm::ResidualResamplingAlgorithm(ResamplingCriterion resamplingCriterion): ResamplingAlgorithm(resamplingCriterion)
{
}

ResidualResamplingAlgorithm* ResidualResamplingAlgorithm::clone() const
{
	return new ResidualResamplingAlgorithm(*this);
}

std::vector<uint> ResidualResamplingAlgorithm::obtainIndexes(uint n,const VectorXd &weights) const
{
	VectorXd residues(weights.size());
	uint *timesToBeResampled = new uint[weights.size()];

	uint nDeterministicParticles = 0;
	for(uint iWeight=0;iWeight<weights.size();iWeight++)
	{
		timesToBeResampled[iWeight] = uint(n*weights(iWeight));
		nDeterministicParticles += timesToBeResampled[iWeight];
		residues(iWeight) = double(n)*weights(iWeight) - double(timesToBeResampled[iWeight]);
	}
	uint nParticlesFromResidues = n - nDeterministicParticles;
	residues /= double(nParticlesFromResidues);

	vector<uint> indexes = StatUtil::discrete_rnd(nParticlesFromResidues,residues);

	indexes.reserve(n);

	uint vectorIndex = nParticlesFromResidues;
	for(uint iWeight=0;iWeight<weights.size();iWeight++)
		for(uint j=0;j<timesToBeResampled[iWeight];j++,vectorIndex++)
			indexes.push_back(iWeight);

	delete[] timesToBeResampled;

	return indexes;
}
