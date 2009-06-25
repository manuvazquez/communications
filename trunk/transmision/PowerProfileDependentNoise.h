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
#ifndef POWERPROFILEDEPENDENTNOISE_H
#define POWERPROFILEDEPENDENTNOISE_H

#include <Noise.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <DelayPowerProfile.h>

class PowerProfileDependentNoise : public Noise
{
protected:
	tMatrix _matrix;
	double _varianceConstant,_stdDev;
public:
    PowerProfileDependentNoise(int nOutputs, int length, const DelayPowerProfile &powerProfile);

	virtual double StdDevAt(int n) const {return _stdDev;}
    virtual tVector operator [ ](int n) const;
    virtual void SetSNR(int SNR, double alphabetVariance);
	virtual void Print() const { cout << _matrix;}

};

#endif
