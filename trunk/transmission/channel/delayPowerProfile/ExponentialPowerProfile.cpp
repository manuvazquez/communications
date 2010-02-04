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

ExponentialPowerProfile::ExponentialPowerProfile(int nOutputs, int nInputs, double T, double threshold): DelayPowerProfile(nOutputs, nInputs)
{
	double power,delay = 0.0;
	double normConst = 0.0;

	while((power=exp(-delay/1e-6))>threshold)
	{
		_tapsPowers.push_back(power);
		delay += T;
		normConst += power;
	}

	for(uint i=0;i<_tapsPowers.size();i++)
		// normalization
		_tapsPowers[i] /= normConst;

	GenerateMatrices();
}

ExponentialPowerProfile::ExponentialPowerProfile(int nOutputs, int nInputs, int m, double tRMS, double T): DelayPowerProfile(nOutputs, nInputs)
{
	double power,delay = 0.0;
	double normConst = 0.0;

	int i;
	for(i=0;i<m;i++)
	{
		power = (1.0/tRMS)*exp(-(1.0/tRMS)*delay);
		_tapsPowers.push_back(power);
		delay += T;
		normConst += power;
	}

	for(i=0;i<m;i++)
		// normalization
		_tapsPowers[i] /= normConst;

	std::vector<double> _amplitudesBak = _tapsPowers;
	for(uint i=0;i<_amplitudesBak.size();i++)
		_tapsPowers[_tapsPowers.size()-1-i] = _amplitudesBak[i];

	GenerateMatrices();
}
