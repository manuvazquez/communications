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
#include "ContinuousPowerProfile.h"

using namespace std;

ContinuousPowerProfile::ContinuousPowerProfile(int nRx, int nTx, vector<double> differentialDelays, vector<double> powers): DelayPowerProfile(nRx, nTx),_continuousDelays(differentialDelays.size()),_continuousPowers(powers.size())
{
	if(differentialDelays.size()!=powers.size())
		throw RuntimeException("ContinuousPowerProfile::ContinuousPowerProfile: numbers of delays and powers differ.");

	double accumulatedDelay = _continuousDelays[0] = differentialDelays[0];
	_continuousPowers[0] = pow(10.0,powers[0]/10.0);
	for(uint i=1;i<differentialDelays.size();i++)
	{
		accumulatedDelay += differentialDelays[i];
		_continuousDelays[i] = accumulatedDelay;
		_continuousPowers[i] = pow(10.0,powers[i]/10.0);
	}
}
