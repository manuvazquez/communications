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
#include "PSPvsSMCSystem.h"

PSPvsSMCSystem::PSPvsSMCSystem()
{

    nSurvivors = 2;

    adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;

    velocity = 50/3.6; // (m/s)
    carrierFrequency = 2e9; // (Hz)
    symbolRate = 500e3; // (Hz)
    T = 1.0/symbolRate; // (s)

    vector<double> differentialDelays,powers;

//     // suburban macro
//     differentialDelays.push_back(0);differentialDelays.push_back(0.1408e-6);differentialDelays.push_back(0.0626e-6);
//     differentialDelays.push_back(0.4015e-6);differentialDelays.push_back(1.3820e-6);differentialDelays.push_back(2.8280e-6);
//     powers.push_back(0);powers.push_back(-2.6682);powers.push_back(-6.2147);
//     powers.push_back(-10.4132);powers.push_back(-16.4735);powers.push_back(-22.1898);

//     if(m==6)
//     {
//         // urban macro
//         differentialDelays.push_back(0);differentialDelays.push_back(0.3600e-6);differentialDelays.push_back(0.2527e-6);
//         differentialDelays.push_back(1.0387e-6);differentialDelays.push_back(2.7300e-6);differentialDelays.push_back(4.5977e-6);
//         powers.push_back(0);powers.push_back(-2.2204);powers.push_back(-1.7184);
//         powers.push_back(-5.1896);powers.push_back(-9.0516);powers.push_back(-12.5013);
//     }else if(m==3)
//     {
//         // urban micro
//         differentialDelays.push_back(0);differentialDelays.push_back(0.2840e-6);differentialDelays.push_back(0.2047e-6);
//         differentialDelays.push_back(0.6623e-6);differentialDelays.push_back(0.8066e-6);differentialDelays.push_back(0.9227e-6);
//         powers.push_back(0);powers.push_back(-1.2661);powers.push_back(-2.7201);
//         powers.push_back(-4.2973);powers.push_back(-6.0140);powers.push_back(-8.4306);
//     }else
//         throw RuntimeException("TesisVariablesAuxiliaresCanalBesselSystem::TesisVariablesAuxiliaresCanalBesselSystem: memory has to be 3 or 6.");

//     powerProfile = new ConstantMeanDSPowerProfile(L,N,differentialDelays,powers,T);
//     powerProfile = new ExponentialPowerProfile(L,N,m,1.8e-6,T);
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

    kalmanEstimator = new KalmanEstimator(powerProfile->means_eigen(),powerProfile->variances_eigen(),N,ARcoefficients,ARvariance);
}

PSPvsSMCSystem::~PSPvsSMCSystem()
{
    delete kalmanEstimator;
    delete powerProfile;
}

void PSPvsSMCSystem::AddAlgorithms()
{
    algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,iLastSymbolVectorToBeDetected+d,ARcoefficients[0],nSurvivors));

    algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen()));
}

void PSPvsSMCSystem::BuildChannel()
{
//  channel = new BesselChannel(N,L,m,symbols.cols(),velocity,carrierFrequency,T,*(dynamic_cast<ContinuousPowerProfile*> (powerProfile)));
    channel = new BesselChannel(N,L,m,symbols.cols(),velocity,carrierFrequency,T,*powerProfile);
}

void PSPvsSMCSystem::BeforeEndingFrame(int iFrame)
{
    SMCSystem::BeforeEndingFrame(iFrame);
    Util::scalarToOctaveFileStream(velocity,"velocity",f);
    Util::scalarToOctaveFileStream(carrierFrequency,"carrierFrequency",f);
    Util::scalarToOctaveFileStream(symbolRate,"symbolRate",f);

    Util::scalarToOctaveFileStream(nSurvivors,"nSurvivors",f);
}
