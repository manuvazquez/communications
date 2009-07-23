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
#include "PowerProfileDependentNoise.h"

PowerProfileDependentNoise::PowerProfileDependentNoise(int nOutputs, int length, const DelayPowerProfile &powerProfile): Noise(nOutputs, length),_matrix(StatUtil::randnMatrix(_nOutputs,_length,0.0,1.0)),_stdDev(1.0)
{
	tMatrix variancesMatrix = powerProfile.variances();
	int i,j;
	double variancesSum = 0.0;
	for(i=0;i<variancesMatrix.rows();i++)
		for(j=0;j<variancesMatrix.cols();j++)
			variancesSum += variancesMatrix(i,j);

	_varianceConstant = variancesSum/double(_nOutputs);
}

tVector PowerProfileDependentNoise::operator [ ](int n) const
{
	tVector res(_nOutputs);
	for(int i=0;i<_nOutputs;i++)
		res(i) = _matrix(i,n);
	return res;
}

void PowerProfileDependentNoise::setSNR(int SNR, double alphabetVariance)
{
	double newStdDev = sqrt(pow(10.0,((double)-SNR)/10.0)*alphabetVariance*_varianceConstant);

	_matrix *= (newStdDev/_stdDev);
	_stdDev = newStdDev;
}
