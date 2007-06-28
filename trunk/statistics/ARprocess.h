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
#ifndef ARPROCESS_H
#define ARPROCESS_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <types.h>
#include <Util.h>
#include <Random.h>
#include <StatUtil.h>
#include <lapackpp/gmd.h>
#include <lapackpp/laslv.h>

class ARprocess{

private:
	vector<double> _coefficients;
	double _noiseVariance;
	double _noiseMean;
	int _nCoefficients, _rows, _columns, _iNextMatrix;
	int _iterationsForConvergence;
	tMatrix **_buffer;

	void CommonConstructorsCode(const tMatrix &seed);

public:
	ARprocess(tMatrix seed,vector<double> coefficients,double noiseVariance);
	ARprocess(tMatrix seed,int order,double velocity,double carrierFrequency,double T);
	ARprocess(const ARprocess &arprocess);
	~ARprocess();

	tMatrix NextMatrix();
	static vector<double> ParametersFromYuleWalker(int order,double velocity,double carrierFrequency,double T,double &noiseVariance);

	int Rows() {return _rows;}
	int Cols() {return _columns;}
};

#endif
