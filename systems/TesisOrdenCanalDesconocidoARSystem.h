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
#ifndef TESISORDENCANALDESCONOCIDOARSYSTEM_H
#define TESISORDENCANALDESCONOCIDOARSYSTEM_H

#include <TesisOrdenCanalDesconocidoSystem.h>

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

class TesisOrdenCanalDesconocidoARSystem : public TesisOrdenCanalDesconocidoSystem
{
protected:
    double channelVariance;

//     virtual void AddAlgorithms();
    virtual void buildChannel();
    virtual void beforeEndingFrame();
public:
    TesisOrdenCanalDesconocidoARSystem();
};

#endif
