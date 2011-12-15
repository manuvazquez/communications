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
#include "TesisOrdenCanalDesconocidoARSystem.h"

TesisOrdenCanalDesconocidoARSystem::TesisOrdenCanalDesconocidoARSystem()
 : TesisOrdenCanalDesconocidoSystem()
{
    channelVariance = 1.0;
//     powerProfile = new FlatPowerProfile(L,N,m,channelVariance);
}

void TesisOrdenCanalDesconocidoARSystem::buildSystemSpecificVariables()
{
    _channel = new ARchannel(_N,_L,_m,_symbols.cols(),ARprocess(_powerProfile->generateChannelMatrix(_randomGenerator),_ARcoefficients,_ARvariance));
}

void TesisOrdenCanalDesconocidoARSystem::saveFrameResults()
{
    TesisOrdenCanalDesconocidoSystem::saveFrameResults();
    Octave::toOctaveFileStream(channelVariance,"channelVariance",_f);
}
