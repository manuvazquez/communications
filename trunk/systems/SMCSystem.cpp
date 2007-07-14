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
#include "SMCSystem.h"

SMCSystem::SMCSystem()
 : BaseSystem(),ARcoefficients(1)
{
    nParticles = 200;
    resamplingRatio = 0.9;

    // back and forward smoothing
    c = 0;
    e = d;

    // system parameters for generating the AR process
    ARprocessOrder = 2;
    velocity = 50; // (Km/h)
    carrierFrequency = 2e9; // (Hz)
    symbolRate = 500e3; // (Hz)
    T = 1.0/symbolRate; // (s)

    // AR process parameters
    ARcoefficients[0] = 0.99999;
    double ARvariance=0.0001;

    vector<double> differentialDelays,powers;
    // suburban macro
//  differentialDelays.push_back(0);differentialDelays.push_back(0.1408e-6);differentialDelays.push_back(0.0626e-6);
//  differentialDelays.push_back(0.4015e-6);differentialDelays.push_back(1.3820e-6);differentialDelays.push_back(2.8280e-6);
//  powers.push_back(0);powers.push_back(-2.6682);powers.push_back(-6.2147);
//  powers.push_back(-10.4132);powers.push_back(-16.4735);powers.push_back(-22.1898);

    // urban macro
//  differentialDelays.push_back(0);differentialDelays.push_back(0.3600e-6);differentialDelays.push_back(0.2527e-6);
//  differentialDelays.push_back(1.0387e-6);differentialDelays.push_back(2.7300e-6);differentialDelays.push_back(4.5977e-6);
//  powers.push_back(0);powers.push_back(-2.2204);powers.push_back(-1.7184);
//  powers.push_back(-5.1896);powers.push_back(-9.0516);powers.push_back(-12.5013);

    // urban micro
    differentialDelays.push_back(0);differentialDelays.push_back(0.2840e-6);differentialDelays.push_back(0.2047e-6);
    differentialDelays.push_back(0.6623e-6);differentialDelays.push_back(0.8066e-6);differentialDelays.push_back(0.9227e-6);
    powers.push_back(0);powers.push_back(-1.2661);powers.push_back(-2.7201);
    powers.push_back(-4.2973);powers.push_back(-6.0140);powers.push_back(-8.4306);

    powerProfile = new ConstantMeanDSPowerProfile(L,N,differentialDelays,powers,T);

    kalmanEstimator = new KalmanEstimator(powerProfile->Means(),powerProfile->Variances(),N,ARcoefficients[0],ARvariance);
    knownSymbolsKalmanEstimator = new KnownSymbolsKalmanEstimator(powerProfile->Means(),powerProfile->Variances(),N,ARcoefficients[0],ARvariance,symbols,preambleLength);

    // ...and linear detectors
    mmseDetectorLarge = new MMSEDetector(L*(c+d+1),N*(m+c+d),alphabet->Variance(),N*(d+1));
    mmseDetectorSmall = new MMSEDetector(L*(c+d+1),N*(d+1),alphabet->Variance(),N*(d+1));
    mmseDetectorXL = new MMSEDetector(L*(c+e+1),N*(e+1),alphabet->Variance(),N*(d+1),0);
    decorrelatorDetector = new DecorrelatorDetector(L*(c+d+1),N*(d+1),alphabet->Variance());

    // always the same resampling criterion and algorithms
    ResamplingCriterion criterioRemuestreo(resamplingRatio);

    algoritmoRemuestreo = new ResidualResamplingAlgorithm(criterioRemuestreo);

    firstSampledChannelMatrixVariance = 0.0;
}


SMCSystem::~SMCSystem()
{
  delete channel;
  delete powerProfile;
  delete kalmanEstimator;
  delete knownSymbolsKalmanEstimator;
  delete mmseDetectorLarge;
  delete mmseDetectorSmall;
  delete mmseDetectorXL;
  delete decorrelatorDetector;
  delete algoritmoRemuestreo;
}

void SMCSystem::BuildChannel()
{
  channel = new BesselChannel(N,L,m,symbols.cols(),velocity,carrierFrequency,T,*powerProfile);

  KnownChannelChannelMatrixEstimator knownChannelEstimator(*channel,preambleLength,N);
}

void SMCSystem::AddAlgorithms()
{
//     algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances()));

//     algorithms.push_back(new SISoptAlgorithm ("SIS opt",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances()));

//     algorithms.push_back(new TriangularizationBasedSMCAlgorithm("Cholesky",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],ARvariance));

    algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF (MMSE)",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,mmseDetectorSmall,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

//     algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF (Decorrelator)",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,decorrelatorDetector,preamble,c,d,d,nParticles,algoritmoRemuestreo,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

//     algorithms.push_back(new ViterbiAlgorithm("Viterbi",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d));

//     algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimator("Kalman Filter (Known Symbols)",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,symbols));

//     algorithms.push_back(new LinearFilterBasedAlgorithm("Kalman Filter + MMSE",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,c,d,mmseDetectorSmall,ARcoefficients[0]));

//     algorithms.push_back(new LinearFilterBasedAlgorithm("Kalman Filter (known symbols) + MMSE",*alphabet,L,N,lastSymbolVectorInstant,m,knownSymbolsKalmanEstimator,preamble,c,d,mmseDetectorSmall,ARcoefficients[0]));
}
