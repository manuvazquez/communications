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
#include "ExponentialPowerProfile.h"

ExponentialPowerProfile::ExponentialPowerProfile(int nRx, int nTx, double T, double threshold): DelayPowerProfile(nRx, nTx)
{
	double power,delay = 0.0;
	double normConst = 0.0;

	while((power=exp(-delay/1e-6))>threshold)
	{
		_amplitudes.push_back(power);
		delay += T;
		normConst += power;
	}

	for(uint i=0;i<_amplitudes.size();i++)
		// normalization
		_amplitudes[i] /= normConst;

	GenerateMatrices();
}
