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
#ifndef TESISORDENCANALDESCONOCIDOSYSTEM_H
#define TESISORDENCANALDESCONOCIDOSYSTEM_H

#include <ChannelOrderEstimationSystem.h>

/**
    @author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <RLSEstimator.h>
#include <RMMSEDetector.h>
#include <WithoutReplacementResamplingAlgorithm.h>
#include <BestParticlesResamplingAlgorithm.h>
#include <MultinomialResamplingAlgorithm.h>
#include <FlatPowerProfile.h>
#include <ExponentialPowerProfile.h>
#include <MLSDmAlgorithm.h>
#include <PSPAlgorithm.h>
#include <TimeInvariantChannel.h>
#include <CMEBasedAlgorithm.h>
#include <TimeVaryingChannelCMEbasedAlgorithm.h>
#include <TransitionCriterion.h>
#include <MaximumProbabilityCriterion.h>
#include <UniformRelatedCriterion.h>
#include <ISIS.h>
#include <USIS.h>
#include <USIS2SISAlgorithm.h>
#include <APPbasedChannelOrderEstimator.h>
#include <LinearFilterBasedSMCAlgorithm.h>

class TesisOrdenCanalDesconocidoSystem : public ChannelOrderEstimationSystem
{
protected:
    int nSurvivors;
    bool adjustParticlesNumberFromSurvivors,adjustSurvivorsFromParticlesNumber;

    double forgettingFactor;
    double forgettingFactorDetector;

    // vectors of channel estimators and linear detectors for unknown channel order algorithms
    vector<ChannelMatrixEstimator *> RLSchannelEstimators;
    vector<ChannelMatrixEstimator *> kalmanChannelEstimators;
    vector<ChannelMatrixEstimator *> noForgetRLSchannelEstimators;
    vector<LinearDetector *> RMMSElinearDetectors;

    RLSEstimator *rlsEstimator;
    RMMSEDetector *rmmseDetector;

    ResamplingAlgorithm *withoutReplacementResamplingAlgorithm,*bestParticlesResamplingAlgorithm,*multinomialResamplingAlgorithm;

    APPbasedChannelOrderEstimator *channelOrderEstimator;
    TransitionCriterion *USISmaximumProbabilityCriterion,*USISuniformRelatedCriterion;

    KalmanEstimator *kalmanEstimator;

    virtual void AddAlgorithms();
    virtual void BeforeEndingFrame(int iFrame);
public:
    TesisOrdenCanalDesconocidoSystem();

    ~TesisOrdenCanalDesconocidoSystem();
};

#endif
