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
#include "ML_SMCAlgorithm.h"

ML_SMCAlgorithm::ML_SMCAlgorithm(string name, Alphabet alphabet, ChannelMatrixEstimator& channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingCriterion resamplingCriterion): SMCAlgorithm(name, alphabet, channelEstimator, preamble, smoothingLag, nParticles, resamplingCriterion)
{
}


void ML_SMCAlgorithm::Process(tMatrix observations, vector< double > noiseVariances)
{
// 	cout << "En Process" << endl;
	int k,iSmoothingVector;
	int iSmoothingLag,iParticle,iSampledVector;
	vector<tSymbol> testedVector(_N),testedSmoothingVector(_N*_d),sampledVector(_N);
	double auxLikelihoodsProd;
	KalmanEstimator *channelEstimatorClone;
	

	// it selects all rows in the symbols Matrix
	tRange allSymbolRows(0,_N-1);

	// it includes all symbol vectors involved in the smoothing
	tMatrix smoothingSymbolVectors(_N,_m+_d);

	int nObservations = observations.cols();
	int nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);
	int nSmoothingVectors = (int) pow((double)_alphabet.Length(),(double)(_N*_d));
	
	// a likelihood is computed for every possible symbol vector
	tVector likelihoods(nSymbolVectors);

	// and from it, a probability
	tVector probabilities(nSymbolVectors);

	// for each time instant
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_endDetectionTime;iObservationToBeProcessed++)
	{
		for(iParticle=0;iParticle<_nParticles;iParticle++)
		{
			// the m-1 already detected symbol vectors are copied into the matrix
			smoothingSymbolVectors(allSymbolRows,*(new tRange(0,_m-1))).inject(_detectedSymbols[iParticle](allSymbolRows,*(new tRange(iObservationToBeProcessed-_m+1,iObservationToBeProcessed))));

			for(int iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.IntToSymbolsArray(iTestedVector,testedVector);
				
				// current tested vector is copied in the m-th position
				for(k=0;k<_N;k++)
					smoothingSymbolVectors(k,_m-1) = testedVector[k];

			} // for(int iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

		} // for(iParticle=0;iParticle<_nParticles;iParticle++)
		
	} // for each time instant
}

