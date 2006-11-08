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
#include "ARprocess.h"

using namespace std;

ARprocess::ARprocess(tMatrix seed,vector<double> coefficients,double noiseVariance):_coefficients(coefficients),_noiseVariance(noiseVariance),_noiseMean(0),_rows(seed.rows()),_columns(seed.cols()),_nCoefficients(coefficients.size()),_buffer(new tMatrix*[_nCoefficients]),_iterationsForConvergence(200)
{
	//the buffer is filled
	int i,j;

	_buffer[0] = new tMatrix(seed);
	for(i=1;i<_nCoefficients;i++)
	{
		_buffer[i] = new tMatrix(_rows,_columns);
		*(_buffer[i]) = 0.0;
		for(j=0;j<i;j++)
			// _buffer[i] = _buffer[i] + _buffer[j]*_coefficients[j];
			Util::Add(*(_buffer[i]),*(_buffer[j]),*(_buffer[i]),1.0,_coefficients[j]);

		tMatrix noise = StatUtil::RandnMatrix(_rows,_columns,_noiseMean,_noiseVariance);

		//_buffer[i] = _buffer[i] + noise;
		Util::Add(*(_buffer[i]),noise,*(_buffer[i]));
	}

	//convergence
	tMatrix aux(_rows,_columns);
	for(i=0;i<_iterationsForConvergence;i++)
	{
		aux = 0.0;
		for(j=0;j<_nCoefficients;j++)
			// aux = aux + _coefficients[j]*_buffer[(i+nCoefficientes-1-j) % nCoefficientes];
			Util::Add(aux,*(_buffer[(i+_nCoefficients-1-j) % _nCoefficients]),aux,1.0,_coefficients[j]);

		tMatrix noise = StatUtil::RandnMatrix(_rows,_columns,_noiseMean,_noiseVariance);

		// _buffer[i % _nCoefficients] = aux + noise;
		Util::Add(aux,noise,*(_buffer[i % _nCoefficients]));
	}

	// the index of the next matrix is gonna be returned is kept
	_iNextMatrix = i;
}

ARprocess::ARprocess(const ARprocess &arprocess):_coefficients(arprocess._coefficients),_noiseVariance(arprocess._noiseVariance),_noiseMean(arprocess._noiseMean),_nCoefficients(arprocess._nCoefficients),_rows(arprocess._rows),_columns(arprocess._columns),_iNextMatrix(arprocess._iNextMatrix),_iterationsForConvergence(arprocess._iterationsForConvergence),_buffer(new tMatrix*[_nCoefficients])/*,_randomGenerator(arprocess._randomGenerator)*/
{
	for(int i=0;i<_nCoefficients;i++)
	{
		_buffer[i] = new tMatrix(_rows,_columns);
		*(_buffer[i]) = *(arprocess._buffer[i]);
	}
}

ARprocess::~ARprocess()
{
	for(int i=0;i<_nCoefficients;i++)
		delete _buffer[i];
	delete[] _buffer;
}

tMatrix ARprocess::NextMatrix()
{
	tMatrix aux(_rows,_columns);

	aux = 0.0;
	for(int j=0;j<_nCoefficients;j++)
		// aux = aux + _coefficients[j]*_buffer[(i+nCoefficientes-1-j) % nCoefficientes];
		Util::Add(aux,*(_buffer[(_iNextMatrix+_nCoefficients-1-j) % _nCoefficients]),aux,1.0,_coefficients[j]);

	tMatrix noise = StatUtil::RandnMatrix(_rows,_columns,_noiseMean,_noiseVariance);

	// _buffer[i % _nCoefficients] = aux + noise;
	Util::Add(aux,noise,*(_buffer[_iNextMatrix % _nCoefficients]));

	return *(_buffer[_iNextMatrix++ % _nCoefficients]);
}


