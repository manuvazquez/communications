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

void TesisComplejidadReducidaBesselNumeroParticulasSystem::addAlgorithms()
{
    for(uint iNparticles=0;iNparticles<particlesNumbers.size();iNparticles++)
    {
        cout << "running algorithms with " << particlesNumbers[iNparticles] << endl;

        // ---------------------------------------------------------- con variables auxiliares ----------------------------------------------------
		std::stringstream choleskyName;
		choleskyName << "Cholesky: " << particlesNumbers[iNparticles] << " particles";
        _algorithms.push_back(new TriangularizationBasedSMCAlgorithm(choleskyName.str(),*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,particlesNumbers[iNparticles],algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_ARcoefficients[0],_ARvariance));

        // aquí restamos la contribución de los símbolos anteriores (el true al final) por lo que se debe usar "mmseDetectorSmall"
		std::stringstream MMSEMKFname;
		MMSEMKFname << "MKF (MMSE): " << particlesNumbers[iNparticles] << " particles";
        _algorithms.push_back(new LinearFilterBasedMKFAlgorithm(MMSEMKFname.str(),*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,mmseDetectorSmall,_preamble,c,_d,_d,particlesNumbers[iNparticles],algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_ARcoefficients[0],firstSampledChannelMatrixVariance,_ARvariance,true));

        // aquí restamos la contribución de los símbolos anteriores (el true al final) por lo que se debe usar "mmseDetectorSmall"
		std::stringstream MMSEDecName;
		MMSEDecName << "MKF (Decorrelator): " << particlesNumbers[iNparticles] << " particles";
        _algorithms.push_back(new LinearFilterBasedMKFAlgorithm(MMSEDecName.str(),*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,decorrelatorDetector,_preamble,c,_d,_d,particlesNumbers[iNparticles],algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_ARcoefficients[0],firstSampledChannelMatrixVariance,_ARvariance,true));

        // ------------------------------------------------ estimacion conjunta del canal y los datos ---------------------------------------------
		std::stringstream RLSDSISname;
		RLSDSISname << "RLS-D-SIS: " << particlesNumbers[iNparticles] << " particles";
        _algorithms.push_back(new LinearFilterBasedSMCAlgorithm(RLSDSISname.str(),*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,rlsEstimator,rmmseDetector,_preamble,c,_d,_d,particlesNumbers[iNparticles],algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_ARcoefficients[0],firstSampledChannelMatrixVariance,_ARvariance));

		std::stringstream LMSDSISname;
		LMSDSISname << "LMS-D-SIS: " << particlesNumbers[iNparticles] << " particles";
        _algorithms.push_back(new LinearFilterBasedSMCAlgorithm(LMSDSISname.str(),*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,lmsEstimator,rmmseDetector,_preamble,c,_d,_d,particlesNumbers[iNparticles],algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_ARcoefficients[0],firstSampledChannelMatrixVariance,_ARvariance));

        // -------------------------------------------------------------- algoritmos comunes ------------------------------------------------------
//         algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,particlesNumbers[iNparticles],algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen()));

//         sprintf(algorithmName,"SIS opt: %d particles)",particlesNumbers[iNparticles]);
//         algorithms.push_back(new SISoptAlgorithm (algorithmName,*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,particlesNumbers[iNparticles],algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen()));
    }

    _algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,_iLastSymbolVectorToBeDetected+_d,nSurvivors));
    _algorithms.push_back(new ViterbiAlgorithm("Viterbi",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (_channel)),_preamble,_d));
}

void TesisComplejidadReducidaBesselNumeroParticulasSystem::saveFrameResults()
{
    TesisComplejidadReducidaBesselSystem::saveFrameResults();
    Octave::toOctaveFileStream(particlesNumbers,"particlesNumbers",_f);
}
