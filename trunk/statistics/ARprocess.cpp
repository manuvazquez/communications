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

#define DEBUG2
#define DEBUG3
// #define DEBUGDOPPLER

using namespace std;

ARprocess::ARprocess(tMatrix seed,vector<double> coefficients,double noiseVariance):_coefficients(coefficients),_noiseVariance(noiseVariance),_noiseMean(0),_nCoefficients(coefficients.size()),_rows(seed.rows()),_columns(seed.cols()),_buffer(new tMatrix*[_nCoefficients])
{
	_buffer[0] = new tMatrix(seed);
	cout << "seed es " << endl << seed;

	CommonConstructorsCode();
}

ARprocess::ARprocess(tMatrix seed,int order,double velocity,double carrierFrequency,double T):_rows(seed.rows()),_columns(seed.cols())
{
	_coefficients = ParametersFromYuleWalker(order,velocity,carrierFrequency,T,_noiseVariance);

#ifdef DEBUG2
	cout << "Los coeficientes son " << endl;
	Util::Print(_coefficients);
	cout << "La varianza es " << _noiseVariance << endl;
#endif

	_noiseMean = 0.0;
	_nCoefficients = _coefficients.size();
	_buffer = new tMatrix*[_nCoefficients];
	_buffer[0] = new tMatrix(seed);

	cout << "seed es " << endl << seed;

	CommonConstructorsCode();
}

void ARprocess::CommonConstructorsCode()
{
	int i,j;

	_iterationsForConvergence = 1000;

	//the buffer is filled
	for(i=1;i<_nCoefficients;i++)
	{
		_buffer[i] = new tMatrix(_rows,_columns);
		*(_buffer[i]) = 0.0;

// 		for(j=0;j<i;j++)
// 			// _buffer[i] = _buffer[i] + _buffer[j]*_coefficients[j];
// 			Util::Add(*(_buffer[i]),*(_buffer[j]),*(_buffer[i]),1.0,_coefficients[j]);
//
// 		tMatrix noise = StatUtil::RandnMatrix(_rows,_columns,_noiseMean,_noiseVariance);
//
// 		//_buffer[i] = _buffer[i] + noise;
// 		Util::Add(*(_buffer[i]),noise,*(_buffer[i]));

		tMatrix noise = StatUtil::RandnMatrix(_rows,_columns,_noiseMean,_noiseVariance);

		//_buffer[i] = _buffer[i] + noise;
		Util::Add(*(_buffer[i-1]),noise,*(_buffer[i]));

#ifdef DEBUG3
		cout << "*(_buffer["<< i <<"])" << endl << *(_buffer[i]) << endl;
#endif
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

#ifdef DEBUG3
		cout << "convergiendo" << endl << *(_buffer[i % _nCoefficients]);
		cout << "Una tecla..."; getchar();
#endif
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

#ifdef DEBUG
	cout << "Otra matriz" << endl << *(_buffer[_iNextMatrix++ % _nCoefficients]);
	cout << "Una tecla..."; getchar();
#endif

	return *(_buffer[_iNextMatrix++ % _nCoefficients]);
}

vector<double> ARprocess::ParametersFromYuleWalker(int order,double velocity,double carrierFrequency,double T,double &noiseVariance)
{
	const double c = 3e8;

	double dopplerFrequency = velocity/(c/carrierFrequency);
#ifdef DEBUG
	cout << "dopplerFrequency = " << dopplerFrequency << endl;
#endif
	double normDopplerFrequency = T*dopplerFrequency;
#ifdef DEBUGDOPPLER
	cout << "normDopplerFrequency = " << normDopplerFrequency << endl;
	normDopplerFrequency = 0.02;
#endif

	tMatrix autocorrelationsMatrix(order,order);
	tVector autocorrelationsVector(order);

	vector<double> autocorrelations(order+1);
	for(int i=0;i<=order;i++)
		autocorrelations[i] = jn(0,2.0*M_PI*normDopplerFrequency*double(i));

#ifdef DEBUG
	cout << "Todas las correlaciones" << endl;
	Util::Print(autocorrelations);
	cout << "jn(0,2.0): " << jn(0,2.0) << endl;
#endif

	for(int m=0;m<order;m++)
	{
		for(int k=0;k<order;k++)
			autocorrelationsMatrix(m,k) = autocorrelations[abs(m-k)];
// 			autocorrelationsMatrix(m,k) = jn(0,2.0*M_PI*normDopplerFrequency*double(abs(m-k)));

// 		autocorrelationsVector(m) = jn(0,2.0*M_PI*normDopplerFrequency*double(m));
		autocorrelationsVector(m) = autocorrelations[m+1];
	}

	tVector coefficients(order);
	LaLinearSolveIP(autocorrelationsMatrix,coefficients,autocorrelationsVector);

	// the variance of the noise will be computed in the loop...
	noiseVariance = autocorrelations[0];

	// ...besides, the lapackpp vector will be copied to a standard c++ vector
	vector<double> coefficientsCppVector(order);
	for(int i=0;i<order;i++)
	{
		noiseVariance -= coefficients(i)*autocorrelations[i+1];
		coefficientsCppVector[i] = coefficients(i);
	}

#ifdef DEBUG
	cout << "noiseVariance is: " << noiseVariance << endl;
	cout << "autocorrelationsMatrix" << endl << autocorrelationsMatrix;
	cout << "autocorrelationsVector" << endl << autocorrelationsVector;
#endif
	return coefficientsCppVector;
}
