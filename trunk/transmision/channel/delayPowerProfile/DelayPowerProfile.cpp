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
#include "DelayPowerProfile.h"

DelayPowerProfile::DelayPowerProfile(int nOutputs,int nInputs):_nOutputs(nOutputs),_nInputs(nInputs),_generatedCoefficientsMean(0.0)
{
}

DelayPowerProfile::~DelayPowerProfile()
{
}

tMatrix DelayPowerProfile::generateChannelMatrix(Random &random)
{
	tMatrix res(_nOutputs,_nInputs*_amplitudes.size());

	for(int i=0;i<res.rows();i++)
		for(int j=0;j<res.cols();j++)
			res(i,j) = random.randn()*sqrt(_variances(i,j)) + _means(i,j);

	return res;
}

void DelayPowerProfile::print() const
{
	for(uint i=0;i<_amplitudes.size();i++)
		std::cout << "amplitude = " << _amplitudes[i] << std::endl;
}

void DelayPowerProfile::GenerateMatrices()
{
	// the memory of the channel is "_amplitudes.size()"
	_means = tMatrix(_nOutputs,_nInputs*_amplitudes.size());
	_means = _generatedCoefficientsMean;

	_variances = tMatrix(_nOutputs,_nInputs*_amplitudes.size());
	tRange rAllRows(0,_nOutputs-1);
	for(uint i=0;i<_amplitudes.size();i++)
		//
		_variances(rAllRows,tRange(i*_nInputs,(i+1)*_nInputs-1)) = _amplitudes[i];
}
