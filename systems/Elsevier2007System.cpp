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
#include "Elsevier2007System.h"

Elsevier2007System::Elsevier2007System()
{
    // the required linear detectors are built
    mmseDetectorLarge = new MMSEDetector(_L*(c+_d+1),_N*(_m+c+_d),_alphabet->variance(),_N*(_d+1));
    mmseDetectorSmall = new MMSEDetector(_L*(c+_d+1),_N*(_d+1),_alphabet->variance(),_N*(_d+1));
//     mmseDetectorXL = new MMSEDetector(_L*(c+e+1),_N*(e+1),_alphabet->variance(),_N*(_d+1),0);
	mmseDetectorXL = new MMSEDetector(_L*(c+_d+1),_N*(_d+1),_alphabet->variance(),_N*(_d+1),0);
    decorrelatorDetector = new DecorrelatorDetector(_L*(c+_d+1),_N*(_d+1),_alphabet->variance());

	kalmanEstimatedChannel = NULL;
}

Elsevier2007System::~Elsevier2007System()
{
	delete mmseDetectorLarge;
	delete mmseDetectorSmall;
	delete mmseDetectorXL;
	delete decorrelatorDetector;
	delete kalmanEstimatedChannel;
}

void Elsevier2007System::addAlgorithms()
{
	delete kalmanEstimatedChannel;
	 kalmanEstimatedChannel = new EstimatedMIMOChannel(_N,_L,_m,_symbols.cols(),_preambleLength,kalmanEstimator,_symbols,_observations,_noise->variances());

//     algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen()));

//     algorithms.push_back(new SISoptAlgorithm ("SIS opt",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,nParticles,algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen()));

//     algorithms.push_back(new TriangularizationBasedSMCAlgorithm("Cholesky",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen(),ARcoefficients[0],ARvariance));

	// si no se resta la contribución de los símbolos anteriores (#define SUBSTRACT_CONTRIBUTION_FROM_KNOWN_SYMBOLS en LinearFilterBasedSMCAlgorithm) entonces
	// se debe usar mmseDetectorLarge
    _algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF (MMSE)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,mmseDetectorLarge,_preamble,c,_d,_d,nParticles,algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_ARcoefficients[0],firstSampledChannelMatrixVariance,_ARvariance));

//     algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF (Decorrelator)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,decorrelatorDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

//     algorithms.push_back(new ViterbiAlgorithm("Viterbi",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d));

// 	algorithms.push_back(new ViterbiAlgorithm("Viterbi (estimated channel)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (kalmanEstimatedChannel)),preamble,d));

//     algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,symbols));

//     algorithms.push_back(new LinearFilterBasedAlgorithm("Kalman Filter + MMSE",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,c,d,mmseDetectorSmall,ARcoefficients[0]));

//     algorithms.push_back(new LinearFilterBasedAlgorithm("Kalman Filter (known symbols) + MMSE",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,knownSymbolsKalmanEstimator,preamble,c,d,mmseDetectorSmall,ARcoefficients[0]));
}
