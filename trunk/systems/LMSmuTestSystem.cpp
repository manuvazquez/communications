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

    velocity = 50/3.6; // (m/s)
    carrierFrequency = 2e9; // (Hz)
    symbolRate = 500e3; // (Hz)
    T = 1.0/symbolRate; // (s)

    vector<double> differentialDelays,powers;

    powerProfile = new FlatPowerProfile(L,N,m,1.0);

    powerProfile->print();

    // check the adjustments for particle and survivor numbers
    if(adjustParticlesNumberFromSurvivors && adjustSurvivorsFromParticlesNumber)
        throw RuntimeException("adjustParticlesNumberFromSurvivors y adjustSurvivorsFromParticlesNumber no pueden ser true a la vez.");

    if(adjustParticlesNumberFromSurvivors)
    {
        nParticles = (int)pow((double)alphabet->length(),N*(m-1))*nSurvivors;
        cout << "Number of particles adjusted to " << nParticles << endl;
    }

    if(adjustSurvivorsFromParticlesNumber)
    {
        cout << "Number of survivors adjusted from " << nSurvivors;
        nSurvivors = int(ceil(double(nParticles)/pow(2.0,double(N*(m-1)))));
        cout << " to " << nSurvivors << endl;
    }

    // estimacion conjunta del canal y los datos
    rmmseDetector = new RMMSEDetector(L*(c+d+1),N*(m+c+d),alphabet->variance(),forgettingFactorDetector,N*(d+1));

    for(uint iMu=0;iMu<musLMS.size();iMu++)
        LMSchannelEstimators.push_back(new NLMSEstimator(powerProfile->means(),N,musLMS[iMu]));
}

LMSmuTestSystem::~LMSmuTestSystem()
{
    delete rmmseDetector;

    for(uint iMu=0;iMu<musLMS.size();iMu++)
        delete LMSchannelEstimators[iMu];

    delete powerProfile;
}

void LMSmuTestSystem::AddAlgorithms()
{
    char algorithmName[ALGORITHM_NAME_MAX_LENGTH];

    for(uint iMu=0;iMu<musLMS.size();iMu++)
    {
        sprintf(algorithmName,"LMS-D-SIS mu = %f",musLMS[iMu]);
        algorithms.push_back(new LinearFilterBasedSMCAlgorithm(algorithmName,*alphabet,L,N,lastSymbolVectorInstant,m,LMSchannelEstimators[iMu],rmmseDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
    }
}

void LMSmuTestSystem::BuildChannel()
{
//  channel = new BesselChannel(N,L,m,symbols.cols(),velocity,carrierFrequency,T,*(dynamic_cast<ContinuousPowerProfile*> (powerProfile)));
    channel = new BesselChannel(N,L,m,symbols.cols(),velocity,carrierFrequency,T,*powerProfile);
}

void LMSmuTestSystem::BeforeEndingFrame(int iFrame)
{
    SMCSystem::BeforeEndingFrame(iFrame);
    Util::ScalarToOctaveFileStream(velocity,"velocity",f);
    Util::ScalarToOctaveFileStream(carrierFrequency,"carrierFrequency",f);
    Util::ScalarToOctaveFileStream(symbolRate,"symbolRate",f);

    Util::ScalarToOctaveFileStream(nSurvivors,"nSurvivors",f);
    Util::ScalarToOctaveFileStream(forgettingFactorDetector,"forgettingFactorDetector",f);
    Util::scalarsVectorToOctaveFileStream(musLMS,"musLMS",f);
}
