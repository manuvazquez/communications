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

    nSurvivors = 1;

    forgettingFactor = 0.99;
    forgettingFactorDetector = 0.95;
    muLMS = 0.01;
    muNLMS = 0.1;

    adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;

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

    powerProfile->Print();

    // check the adjustments for particle and survivor numbers
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

    // variables auxiliares
//     mmseDetectorLarge = new MMSEDetector(L*(c+d+1),N*(m+c+d),alphabet->Variance(),N*(d+1));
    mmseDetectorSmall = new MMSEDetector(L*(c+d+1),N*(d+1),alphabet->Variance(),N*(d+1));
    decorrelatorDetector = new DecorrelatorDetector(L*(c+d+1),N*(d+1),alphabet->Variance());

    // estimacion conjunta del canal y los datos
    rmmseDetector = new RMMSEDetector(L*(c+d+1),N*(m+c+d),alphabet->Variance(),forgettingFactorDetector,N*(d+1));
    rlsEstimator = new RLSEstimator(powerProfile->Means(),N,forgettingFactor);
    lmsEstimator = new LMSEstimator(powerProfile->Means(),N,muLMS);
    nlmsEstimator = new NLMSEstimator(powerProfile->Means(),N,muNLMS);

    kalmanEstimator = new KalmanEstimator(powerProfile->Means(),powerProfile->Variances(),N,ARcoefficients[0],ARvariance);
    knownSymbolsKalmanEstimator = new KnownSymbolsKalmanEstimator(powerProfile->Means(),powerProfile->Variances(),N,ARcoefficients[0],ARvariance,symbols,preambleLength);

    kalmanEstimatedChannel = NULL;
}

TesisComplejidadReducidaSystem::~TesisComplejidadReducidaSystem()
{
//     delete mmseDetectorLarge;
    delete mmseDetectorSmall;
    delete decorrelatorDetector;

    delete rmmseDetector;
    delete rlsEstimator;
    delete lmsEstimator;

    delete kalmanEstimatedChannel;
    delete kalmanEstimator;
    delete knownSymbolsKalmanEstimator;
    delete powerProfile;
}

void TesisComplejidadReducidaSystem::AddAlgorithms()
{
    // ---------------------------------------------------------- con variables auxiliares ----------------------------------------------------
    algorithms.push_back(new TriangularizationBasedSMCAlgorithm("Cholesky",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],ARvariance));

    // aquí restamos la contribución de los símbolos anteriores (el true al final) por lo que se debe usar "mmseDetectorSmall"
    algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF (MMSE)",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,mmseDetectorSmall,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,true));

    // aquí restamos la contribución de los símbolos anteriores (el true al final) por lo que se debe usar "mmseDetectorSmall"
    algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF (Decorrelator)",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,decorrelatorDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,true));

    algorithms.push_back(new LinearFilterBasedAlgorithm("Kalman Filter + MMSE",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,c,d,mmseDetectorSmall,ARcoefficients[0],true));

    algorithms.push_back(new LinearFilterBasedAlgorithm("Kalman Filter (known symbols) + MMSE",*alphabet,L,N,lastSymbolVectorInstant,m,knownSymbolsKalmanEstimator,preamble,c,d,mmseDetectorSmall,ARcoefficients[0],true));

    algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,d,lastSymbolVectorInstant+d,ARcoefficients[0],nSurvivors));

    // ------------------------------------------------ estimacion conjunta del canal y los datos ---------------------------------------------
    algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS",*alphabet,L,N,lastSymbolVectorInstant,m,rlsEstimator,rmmseDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
    algorithms.push_back(new LinearFilterBasedSMCAlgorithm("LMS-D-SIS",*alphabet,L,N,lastSymbolVectorInstant,m,lmsEstimator,rmmseDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
    algorithms.push_back(new LinearFilterBasedSMCAlgorithm("NLMS-D-SIS",*alphabet,L,N,lastSymbolVectorInstant,m,nlmsEstimator,rmmseDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    // -------------------------------------------------------------- algoritmos comunes ------------------------------------------------------
    // common to all simulation algorithms
    delete kalmanEstimatedChannel;
     kalmanEstimatedChannel = new EstimatedMIMOChannel(N,L,m,symbols.cols(),preambleLength,kalmanEstimator,symbols,observaciones,ruido->Variances());

//     algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances()));

    algorithms.push_back(new SISoptAlgorithm ("SIS opt",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances()));

    algorithms.push_back(new ViterbiAlgorithm("Viterbi",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d));

    algorithms.push_back(new ViterbiAlgorithm("Viterbi (estimated channel)",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (kalmanEstimatedChannel)),preamble,d));

    algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,symbols));
}

void TesisComplejidadReducidaSystem::BeforeEndingFrame(int iFrame)
{
    SMCSystem::BeforeEndingFrame(iFrame);

    Util::ScalarToOctaveFileStream(nSurvivors,"nSurvivors",f);
    Util::ScalarToOctaveFileStream(forgettingFactor,"forgettingFactor",f);
    Util::ScalarToOctaveFileStream(forgettingFactorDetector,"forgettingFactorDetector",f);
    Util::ScalarToOctaveFileStream(muLMS,"muLMS",f);
    Util::ScalarToOctaveFileStream(muNLMS,"muNLMS",f);
}
