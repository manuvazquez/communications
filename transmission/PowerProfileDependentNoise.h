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
	MatrixXd _matrix;
	
	/**
	 * @brief it stores those computations involved in the calculus of the variance/std that don't depend on the SNR
	 **/
	double _powerProfileDependentVarianceFactor;
	
	double _stdDev;
	double _alphabetVariance;
	
	virtual double computeStd(const int &SNR) const
	{
		return sqrt(pow(10.0,((double)-SNR)/10.0)*_powerProfileDependentVarianceFactor);
	}
public:
    PowerProfileDependentNoise(double alphabetVariance, uint nOutputs, uint length, const DelayPowerProfile &powerProfile);

	virtual double stdDevAt(uint n) const {return _stdDev;}
    virtual VectorXd at(uint n) const;
    virtual void setSNR(int SNR);
	virtual void print() const { cout << _matrix;}

	static std::string getXMLname() { return "PowerProfileDependentNoise"; }
};

#endif
