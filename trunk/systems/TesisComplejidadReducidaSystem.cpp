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
#include "TesisComplejidadReducidaSystem.h"

TesisComplejidadReducidaSystem::TesisComplejidadReducidaSystem()
{

    nSurvivors = 2;
    adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;

    if(adjustParticlesNumberFromSurvivors && adjustSurvivorsFromParticlesNumber)
        throw RuntimeException("adjustParticlesNumberFromSurvivors y adjustSurvivorsFromParticlesNumber no pueden ser true a la vez.");

    if(adjustParticlesNumberFromSurvivors)
    {
        nParticles = (int)pow((double)alphabet->Length(),N*(m-1))*nSurvivors;
        cout << "Number of particles adjusted to " << nParticles << endl;
    }

    if(adjustSurvivorsFromParticlesNumber)
    {
        cout << "Number of survivors adjusted from " << nSurvivors;
        nSurvivors = int(ceil(double(nParticles)/pow(2.0,double(N*(m-1)))));
        cout << " to " << nSurvivors << endl;
    }

    // the required linear detectors are built
//     mmseDetectorLarge = new MMSEDetector(L*(c+d+1),N*(m+c+d),alphabet->Variance(),N*(d+1));
//     mmseDetectorSmall = new MMSEDetector(L*(c+d+1),N*(d+1),alphabet->Variance(),N*(d+1));
//     mmseDetectorXL = new MMSEDetector(L*(c+e+1),N*(e+1),alphabet->Variance(),N*(d+1),0);
//     decorrelatorDetector = new DecorrelatorDetector(L*(c+d+1),N*(d+1),alphabet->Variance());

    kalmanEstimatedChannel = NULL;
}

TesisComplejidadReducidaSystem::~TesisComplejidadReducidaSystem()
{
//     delete mmseDetectorLarge;
//     delete mmseDetectorSmall;
//     delete mmseDetectorXL;
//     delete decorrelatorDetector;
    delete kalmanEstimatedChannel;
}

void TesisComplejidadReducidaSystem::AddAlgorithms()
{
    delete kalmanEstimatedChannel;
     kalmanEstimatedChannel = new EstimatedMIMOChannel(N,L,m,symbols.cols(),preambleLength,kalmanEstimator,symbols,observaciones,ruido->Variances());

    algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances()));

    algorithms.push_back(new SISoptAlgorithm ("SIS opt",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances()));

//     algorithms.push_back(new TriangularizationBasedSMCAlgorithm("Cholesky",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],ARvariance));

    // si no se resta la contribución de los símbolos anteriores (#define SUBSTRACT_CONTRIBUTION_FROM_KNOWN_SYMBOLS en LinearFilterBasedSMCAlgorithm) entonces
    // se debe usar mmseDetectorLarge
//     algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF (MMSE)",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,mmseDetectorLarge,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

//     algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF (Decorrelator)",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,decorrelatorDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    algorithms.push_back(new ViterbiAlgorithm("Viterbi",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d));

    algorithms.push_back(new ViterbiAlgorithm("Viterbi (estimated channel)",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (kalmanEstimatedChannel)),preamble,d));

    algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,symbols));

//     algorithms.push_back(new LinearFilterBasedAlgorithm("Kalman Filter + MMSE",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,c,d,mmseDetectorSmall,ARcoefficients[0]));
//
//     algorithms.push_back(new LinearFilterBasedAlgorithm("Kalman Filter (known symbols) + MMSE",*alphabet,L,N,lastSymbolVectorInstant,m,knownSymbolsKalmanEstimator,preamble,c,d,mmseDetectorSmall,ARcoefficients[0]));
}
