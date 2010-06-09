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
#ifndef PSPVSPSPBASEDSMCSYSTEM_H
#define PSPVSPSPBASEDSMCSYSTEM_H

#include <SMCSystem.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <KalmanEstimator.h>
#include <ViterbiAlgorithm.h>
#include <PSPAlgorithm.h>
#include <PSPBasedSMCAlgorithm.h>
#include <FlatPowerProfile.h>
#include <WithoutReplacementResamplingAlgorithm.h>
#include <BestParticlesResamplingAlgorithm.h>

class PSPvsPSPBasedSMCSystem : public SMCSystem
{
public:
    PSPvsPSPBasedSMCSystem();

    ~PSPvsPSPBasedSMCSystem();

protected:
    int nSurvivors;
    bool adjustParticlesNumberFromSurvivors;

    KalmanEstimator *kalmanEstimator;
    ResamplingAlgorithm *withoutReplacementResamplingAlgorithm,*bestParticlesResamplingAlgorithm;

    virtual void beforeEndingFrame();
    virtual void addAlgorithms();
    virtual void buildChannel();
};

#endif
