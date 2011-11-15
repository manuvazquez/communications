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
#ifndef NOISE_H
#define NOISE_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <types.h>
#include <StatUtil.h>

class Noise{
protected:
	uint _nOutputs,_length;
public:
    Noise(uint nOutputs,uint length);
	virtual ~Noise() {};

	uint length() const { return _length;}
	uint nOutputs() const { return _nOutputs;}
	virtual void print() const = 0;
	virtual double stdDevAt(int n) const = 0;
    
    virtual VectorXd at(uint n) const = 0;
    
    virtual MatrixXd range(int start,int end) const {throw RuntimeException("Noise::range: not implemented.");}   
    
    
	double VarianceAt(int n) const { double stdDev = stdDevAt(n); return stdDev*stdDev;};
	vector<double> variances() const;
	virtual void setSNR(int SNR) = 0;
};

#endif
