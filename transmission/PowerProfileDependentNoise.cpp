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

PowerProfileDependentNoise::PowerProfileDependentNoise(double alphabetVariance, uint nOutputs, uint length, const DelayPowerProfile &powerProfile): Noise(nOutputs, length),_matrix(StatUtil::randnMatrix(_nOutputs,_length,0.0,1.0)),_stdDev(1.0),_alphabetVariance(alphabetVariance)
{	
// 	_powerProfileDependentVarianceFactor = _alphabetVariance * ( (powerProfile.variances().array() + powerProfile.means().array()*powerProfile.means().array()).sum()/double(_nOutputs) );
	_powerProfileDependentVarianceFactor = _alphabetVariance * ( (powerProfile.variances().array() + powerProfile.means().array()*powerProfile.means().array()).sum()/double(powerProfile.variances().rows()) );
}

VectorXd PowerProfileDependentNoise::at(uint n) const
{
    return _matrix.col(n);
}

void PowerProfileDependentNoise::setSNR(int SNR)
{
	double newStdDev = computeStd(SNR);
	
	_matrix *= (newStdDev/_stdDev);
	_stdDev = newStdDev;
}
