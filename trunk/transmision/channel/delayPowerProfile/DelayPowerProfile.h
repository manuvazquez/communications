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
#ifndef DELAYPOWERPROFILE_H
#define DELAYPOWERPROFILE_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <vector>
#include <Random.h>

class DelayPowerProfile{
protected:
	int _nRx,_nTx;
	std::vector<double> _delays,_amplitudes;
	double _generatedCoefficientsMean;
	tMatrix _means,_variances;

	void GenerateMatrices();
public:
    DelayPowerProfile(int nRx,int nTx);

	virtual ~DelayPowerProfile();

	virtual tMatrix GenerateChannelMatrix(Random &random);
	virtual void Print();
	tMatrix Means() { return _means;}
	tMatrix Variances() {return _variances;}
	int Memory() {return _amplitudes.size();}
};

#endif
