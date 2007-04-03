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

ResidualResamplingAlgorithm* ResidualResamplingAlgorithm::Clone() const
{
	return new ResidualResamplingAlgorithm(*this);
}

std::vector<int> ResidualResamplingAlgorithm::ObtainIndexes(int n,const tVector &weights) const
{
	tVector residues(weights.size());
	int *timesToBeResampled = new int[weights.size()];

	int nDeterministicParticles = 0;
	for(int iWeight=0;iWeight<weights.size();iWeight++)
	{
		timesToBeResampled[iWeight] = n*weights(iWeight);
		nDeterministicParticles += timesToBeResampled[iWeight];
		residues(iWeight) = double(n)*weights(iWeight) - double(timesToBeResampled[iWeight]);
	}
	int nParticlesFromResidues = n - nDeterministicParticles;
	residues *= 1.0/double(nParticlesFromResidues);

	vector<int> indexes = StatUtil::Discrete_rnd(nParticlesFromResidues,residues);

	indexes.reserve(n);

	int vectorIndex = nParticlesFromResidues;
	for(int iWeight=0;iWeight<weights.size();iWeight++)
		for(int j=0;j<timesToBeResampled[iWeight];j++,vectorIndex++)
			indexes.push_back(iWeight);

	delete[] timesToBeResampled;

	return indexes;
}

// void ResidualResamplingAlgorithm::Resample(ParticleFilter* particleFilter, const tVector& weights)
// {
// 	tVector residues(particleFilter->Nparticles());
// 	int *timesToBeResampled = new int[particleFilter->Nparticles()];
//
// 	int nDeterministicParticles = 0;
// 	for(int iParticle=0;iParticle<particleFilter->Nparticles();iParticle++)
// 	{
// 		timesToBeResampled[iParticle] = particleFilter->Nparticles()*weights(iParticle);
// 		nDeterministicParticles += timesToBeResampled[iParticle];
// 		residues(iParticle) = double(particleFilter->Nparticles())*weights(iParticle) - double(timesToBeResampled[iParticle]);
// 		#ifdef DEBUG
// 			cout << "timesToBeResampled[iParticle] = " << timesToBeResampled[iParticle] << endl;
// 			cout << "El residuo es " << residues(iParticle) << endl;
// 		#endif
// 	}
// 	int nParticlesFromResidues = particleFilter->Nparticles() - nDeterministicParticles;
// 	residues *= 1.0/double(nParticlesFromResidues);
//
// 	vector<int> indexes = StatUtil::Discrete_rnd(nParticlesFromResidues,residues);
//
// 	indexes.reserve(particleFilter->Nparticles());
//
// 	int vectorIndex = nParticlesFromResidues;
//
// 	for(int iParticle=0;iParticle<particleFilter->Nparticles();iParticle++)
// 		for(int j=0;j<timesToBeResampled[iParticle];j++,vectorIndex++)
// 			indexes.push_back(iParticle);
//
// 	particleFilter->KeepParticles(indexes);
//
// 	#ifdef DEBUG
// 		cout << "Los pesos son" << endl << weights;
// 		cout << "nº de veces sin residuo" << endl;
// 		Util::Print(timesToBeResampled,particleFilter->Nparticles());
// 		cout << "los residuos" << endl << residues;
// 		cout << "los indices" << endl;
// 		Util::Print(indexes);
// 		cout << "nParticlesFromResidues = " << nParticlesFromResidues << endl;
// 		cout << "Una tecla..."; getchar();
// 	#endif
//
// 	delete[] timesToBeResampled;
// }
