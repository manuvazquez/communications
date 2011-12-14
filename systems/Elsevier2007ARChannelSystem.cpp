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
    _powerProfile = new FlatPowerProfile(_L,_N,_m,channelVariance);

    kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,ARcoefficients,ARvariance);
    knownSymbolsKalmanEstimator = new KnownSymbolsKalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,ARcoefficients,ARvariance,_symbols,_preambleLength);
}


Elsevier2007ARChannelSystem::~Elsevier2007ARChannelSystem()
{
  delete _powerProfile;
  delete kalmanEstimator;
  delete knownSymbolsKalmanEstimator;
}

void Elsevier2007ARChannelSystem::buildSystemSpecificVariables()
{
    _channel = new ARchannel(_N,_L,_m,_symbols.cols(),ARprocess(_powerProfile->generateChannelMatrix(_randomGenerator),ARcoefficients,ARvariance));
}

void Elsevier2007ARChannelSystem::beforeEndingFrame()
{
    Elsevier2007System::beforeEndingFrame();
    Util::scalarToOctaveFileStream(channelVariance,"channelVariance",_f);
}

void Elsevier2007ARChannelSystem::saveFrameResults()
{
    Elsevier2007System::saveFrameResults();
    Octave::toOctaveFileStream(channelVariance,"channelVariance",_f);
}
