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
#ifndef SMCSYSTEM_H
#define SMCSYSTEM_H

#include <BaseSystem.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <DelayPowerProfile.h>
#include <ConstantMeanDSPowerProfile.h>
#include <BesselChannel.h>
#include <KalmanEstimator.h>
#include <KnownChannelChannelMatrixEstimator.h>
#include <KnownSymbolsKalmanEstimator.h>
#include <MMSEDetector.h>
#include <DecorrelatorDetector.h>

#include <DSISoptAlgorithm.h>
#include <LinearFilterBasedSMCAlgorithm.h>
#include <LinearFilterBasedMKFAlgorithm.h>
#include <ViterbiAlgorithm.h>
#include <KnownSymbolsKalmanBasedChannelEstimator.h>
#include <TriangularizationBasedSMCAlgorithm.h>
#include <LinearFilterBasedAlgorithm.h>
#include <SISoptAlgorithm.h>

#include <ResamplingCriterion.h>
#include <ResamplingAlgorithm.h>
#include <ResidualResamplingAlgorithm.h>

class SMCSystem : public BaseSystem
{
protected:
    int nParticles;
    double resamplingRatio;

    // back and forward smoothing
    int c,e;

    // system parameters for generating the AR process
    int ARprocessOrder;
//     double velocity; // (Km/h)
//     double carrierFrequency; // (Hz)
//     double symbolRate; // (Hz)
//     double T; // (s)

//     DelayPowerProfile *powerProfile;

    std::vector<double> ARcoefficients;
    double ARvariance;

//     KalmanEstimator *kalmanEstimator;
//     KnownSymbolsKalmanEstimator *knownSymbolsKalmanEstimator;
//     MMSEDetector *mmseDetectorLarge,*mmseDetectorSmall,*mmseDetectorXL;
//     DecorrelatorDetector *decorrelatorDetector;

    ResamplingAlgorithm *algoritmoRemuestreo;

    double firstSampledChannelMatrixVariance;

    virtual void BeforeEndingFrame(int iFrame);
    virtual void BuildChannel();
//     virtual void AddAlgorithms();
public:
    SMCSystem();

    ~SMCSystem();
};

#endif
