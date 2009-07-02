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
#include "Elsevier2007ARChannelSystem.h"

Elsevier2007ARChannelSystem::Elsevier2007ARChannelSystem()
 : Elsevier2007System()
{
    channelVariance = 1.0;
    powerProfile = new FlatPowerProfile(L,N,m,channelVariance);

    kalmanEstimator = new KalmanEstimator(powerProfile->means(),powerProfile->variances(),N,ARcoefficients,ARvariance);
    knownSymbolsKalmanEstimator = new KnownSymbolsKalmanEstimator(powerProfile->means(),powerProfile->variances(),N,ARcoefficients,ARvariance,symbols,preambleLength);
}


Elsevier2007ARChannelSystem::~Elsevier2007ARChannelSystem()
{
  delete powerProfile;
//   delete channel;
  delete kalmanEstimator;
  delete knownSymbolsKalmanEstimator;
}

void Elsevier2007ARChannelSystem::BuildChannel()
{
    channel = new ARchannel(N,L,m,symbols.cols(),ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
}

void Elsevier2007ARChannelSystem::BeforeEndingFrame(int iFrame)
{
    Elsevier2007System::BeforeEndingFrame(iFrame);
    Util::ScalarToOctaveFileStream(channelVariance,"channelVariance",f);
}
