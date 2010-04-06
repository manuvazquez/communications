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
#include "LinearFilterBasedMKFAlgorithm.h"

LinearFilterBasedMKFAlgorithm::LinearFilterBasedMKFAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, int m, KalmanEstimator* channelEstimator, LinearDetector* linearDetector, MatrixXd preamble, int backwardsSmoothingLag, int smoothingLag, int forwardSmoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const MatrixXd& channelMatrixMean, const MatrixXd& channelMatrixVariances, double ARcoefficient, double samplingVariance, double ARprocessVariance, bool substractContributionFromKnownSymbols): LinearFilterBasedSMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, linearDetector, preamble, backwardsSmoothingLag, smoothingLag, forwardSmoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances, ARcoefficient, samplingVariance, ARprocessVariance,substractContributionFromKnownSymbols)
{
//     _randomParticlesInitilization = true;
}

