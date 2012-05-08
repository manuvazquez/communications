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
#include "LMSmuTestSystem.h"

LMSmuTestSystem::LMSmuTestSystem()
{

    nSurvivors = 2;

    forgettingFactorDetector = 0.95;

//     musLMS.push_back(0.001);musLMS.push_back(0.01);musLMS.push_back(0.1);musLMS.push_back(0.2);musLMS.push_back(0.5);musLMS.push_back(0.9);
    musLMS.push_back(0.8);musLMS.push_back(0.9);musLMS.push_back(0.95);musLMS.push_back(0.99);

    adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;

    vector<double> differentialDelays,powers;

    _powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);

    _powerProfile->print();
	
	adjustParticlesSurvivors(_nParticles,nSurvivors,adjustParticlesNumberFromSurvivors,adjustSurvivorsFromParticlesNumber);

    // estimacion conjunta del canal y los datos
    rmmseDetector = new RMMSEDetector(_L*(_d+1),_N*(_m+_d),_alphabet->variance(),forgettingFactorDetector,_N*(_d+1));

    for(uint iMu=0;iMu<musLMS.size();iMu++)
        LMSchannelEstimators.push_back(new NLMSEstimator(_powerProfile->means(),_N,musLMS[iMu]));
}

LMSmuTestSystem::~LMSmuTestSystem()
{
    delete rmmseDetector;

    for(uint iMu=0;iMu<musLMS.size();iMu++)
        delete LMSchannelEstimators[iMu];

    delete _powerProfile;
}

void LMSmuTestSystem::addAlgorithms()
{
    for(uint iMu=0;iMu<musLMS.size();iMu++)
    {
		std::stringstream algorithmName;
		algorithmName << "LMS-D-SIS mu = " << musLMS[iMu];

        _algorithms.push_back(new LinearFilterBasedSMCAlgorithm(algorithmName.str(),*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,LMSchannelEstimators[iMu],rmmseDetector,_preamble,_d,_nParticles,_resamplingAlgorithm,_powerProfile->means(),_powerProfile->variances(),_ARcoefficients[0],_firstSampledChannelMatrixVariance,_ARvariance));
    }
}

void LMSmuTestSystem::saveFrameResults()
{
    SMCSystem::saveFrameResults();

    Octave::toOctaveFileStream(nSurvivors,"nSurvivors",_f);
    Octave::toOctaveFileStream(forgettingFactorDetector,"forgettingFactorDetector",_f);
    Octave::toOctaveFileStream(musLMS,"musLMS",_f);
}