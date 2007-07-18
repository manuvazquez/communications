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
#ifndef ELSEVIER2007BESSELCHANNELSYSTEM_H
#define ELSEVIER2007BESSELCHANNELSYSTEM_H

#include <Elsevier2007System.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <ExponentialPowerProfile.h>
#include <FlatPowerProfile.h>

class Elsevier2007BesselChannelSystem : public Elsevier2007System
{
protected:
    double velocity; // (Km/h)
    double carrierFrequency; // (Hz)
    double symbolRate; // (Hz)
    double T; // (s)

    virtual void BuildChannel();
    virtual void BeforeEndingFrame(int iFrame);
public:
    Elsevier2007BesselChannelSystem();

    ~Elsevier2007BesselChannelSystem();

};

#endif
