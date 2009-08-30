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
#include "TesisComplejidadReducidaBesselNumeroParticulasSystem.h"

TesisComplejidadReducidaBesselNumeroParticulasSystem::TesisComplejidadReducidaBesselNumeroParticulasSystem(): TesisComplejidadReducidaBesselSystem()
{
    particlesNumbers.push_back(30);particlesNumbers.push_back(50);particlesNumbers.push_back(100);particlesNumbers.push_back(500);particlesNumbers.push_back(1000);
}

void TesisComplejidadReducidaBesselNumeroParticulasSystem::BeforeEndingFrame(int iFrame)
{
    TesisComplejidadReducidaBesselSystem::BeforeEndingFrame(iFrame);
    Util::scalarsVectorToOctaveFileStream(particlesNumbers,"particlesNumbers",f);
}

void TesisComplejidadReducidaBesselNumeroParticulasSystem::AddAlgorithms()
{
    char algorithmName[ALGORITHM_NAME_MAX_LENGTH];

    for(uint iNparticles=0;iNparticles<particlesNumbers.size();iNparticles++)
    {
        cout << "running algorithms with " << particlesNumbers[iNparticles] << endl;

        // ---------------------------------------------------------- con variables auxiliares ----------------------------------------------------
        sprintf(algorithmName,"Cholesky: %d particles)",particlesNumbers[iNparticles]);
        algorithms.push_back(new TriangularizationBasedSMCAlgorithm(algorithmName,*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,particlesNumbers[iNparticles],algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen(),ARcoefficients[0],ARvariance));

        // aquí restamos la contribución de los símbolos anteriores (el true al final) por lo que se debe usar "mmseDetectorSmall"
        sprintf(algorithmName,"MKF (MMSE): %d particles)",particlesNumbers[iNparticles]);
        algorithms.push_back(new LinearFilterBasedMKFAlgorithm(algorithmName,*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,mmseDetectorSmall,preamble,c,d,d,particlesNumbers[iNparticles],algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,true));

        // aquí restamos la contribución de los símbolos anteriores (el true al final) por lo que se debe usar "mmseDetectorSmall"
        sprintf(algorithmName,"MKF (Decorrelator): %d particles)",particlesNumbers[iNparticles]);
        algorithms.push_back(new LinearFilterBasedMKFAlgorithm(algorithmName,*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,decorrelatorDetector,preamble,c,d,d,particlesNumbers[iNparticles],algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,true));

        // ------------------------------------------------ estimacion conjunta del canal y los datos ---------------------------------------------
        sprintf(algorithmName,"RLS-D-SIS: %d particles)",particlesNumbers[iNparticles]);
        algorithms.push_back(new LinearFilterBasedSMCAlgorithm(algorithmName,*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,rlsEstimator,rmmseDetector,preamble,c,d,d,particlesNumbers[iNparticles],algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

        sprintf(algorithmName,"LMS-D-SIS: %d particles)",particlesNumbers[iNparticles]);
        algorithms.push_back(new LinearFilterBasedSMCAlgorithm(algorithmName,*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,lmsEstimator,rmmseDetector,preamble,c,d,d,particlesNumbers[iNparticles],algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

        // -------------------------------------------------------------- algoritmos comunes ------------------------------------------------------
//         algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,particlesNumbers[iNparticles],algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen()));

//         sprintf(algorithmName,"SIS opt: %d particles)",particlesNumbers[iNparticles]);
//         algorithms.push_back(new SISoptAlgorithm (algorithmName,*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,particlesNumbers[iNparticles],algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen()));
    }

    algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,iLastSymbolVectorToBeDetected+d,ARcoefficients[0],nSurvivors));
    algorithms.push_back(new ViterbiAlgorithm("Viterbi",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d));
}
