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
#include "TesisOrdenCanalDesconocidoBesselSystem.h"

TesisOrdenCanalDesconocidoBesselSystem::TesisOrdenCanalDesconocidoBesselSystem()
 : TesisOrdenCanalDesconocidoSystem()
{
    velocity = 50/3.6; // (m/s)
    carrierFrequency = 2e9; // (Hz)
    symbolRate = 500e3; // (Hz)
    T = 1.0/symbolRate; // (s)
}

void TesisOrdenCanalDesconocidoBesselSystem::buildSystemSpecificVariables()
{
    _channel = new BesselChannel(_N,_L,_m,_symbols.cols(),velocity,carrierFrequency,T,*_powerProfile);
	
	// the noise is built here...and it might depend on the channel
	TesisOrdenCanalDesconocidoSystem::buildSystemSpecificVariables();
}

void TesisOrdenCanalDesconocidoBesselSystem::beforeEndingFrame()
{
    TesisOrdenCanalDesconocidoSystem::beforeEndingFrame();
    Util::scalarToOctaveFileStream(velocity,"velocity",_f);
    Util::scalarToOctaveFileStream(carrierFrequency,"carrierFrequency",_f);
    Util::scalarToOctaveFileStream(symbolRate,"symbolRate",_f);
}

