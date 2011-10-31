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
#include "ConstantMeanDSPowerProfile.h"

// #define DEBUG

ConstantMeanDSPowerProfile::ConstantMeanDSPowerProfile(uint nOutputs, uint nInputs, std::vector< double > differentialDelays, std::vector< double > powers, double T): ContinuousPowerProfile(nOutputs, nInputs, differentialDelays, powers)
{
	double quotient;
	int k;
	int nDelays = int(ceil(_continuousDelays[_continuousDelays.size()-1]/T))+1;
	_tapsPowers.resize(nDelays,0.0);

	for(uint i=0;i<_continuousDelays.size();i++)
	{
		quotient = _continuousDelays[i] / T;
		if(quotient == floor(quotient))
		{
			_tapsPowers[int(quotient)] += _continuousPowers[i];
			continue;
		}
		k = int(_continuousDelays[i] / T);
		_tapsPowers[k] += ((k+1)- _continuousDelays[i]/T)*_continuousPowers[i];
		_tapsPowers[k+1] += (_continuousDelays[i]/T - k)*_continuousPowers[i];
	}

	vector<double> _amplitudesBak = _tapsPowers;
	for(uint i=0;i<_amplitudesBak.size();i++)
		_tapsPowers[_tapsPowers.size()-1-i] = _amplitudesBak[i];

	GenerateMatrices();
}
