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

// ARprocess::ARprocess()
// {
// }

ARprocess::ARprocess(tMatrix seed,vector<double> coefficients,double noiseVariance)
{
	this->coefficients = coefficients;
	this->noiseVariance = noiseVariance;
	this->noiseMean = 0;
	this->rows = seed.rows();
	this->columns = seed.cols();
	this->nCoefficients = coefficients.size();
	this->buffer = new tMatrix[nCoefficients];
	this->randomGenerator = new Random(46576);
	this->iterationsForConvergence = 200;

	//the buffer is filled
	int i,j;
	tMatrix noise;
	buffer[0] = seed;
	for(i=1;i<nCoefficients;i++)
	{
		buffer[i] = *(new tMatrix(rows,columns));
		buffer[i] = 0;
		for(j=0;j<i;j++)
			// buffer[i] = buffer[i] + buffer[j]*coefficients[j];
			Util::Add(buffer[i],buffer[j],buffer[i],1.0,coefficients[j]);

// 		tMatrix noise(randomGenerator->randnArray(rows*columns,noiseMean,noiseVariance),rows,columns);
		noise = *(new tMatrix(randomGenerator->randnArray(rows*columns,noiseMean,noiseVariance),rows,columns));

		//buffer[i] = buffer[i] + noise;
		Util::Add(buffer[i],noise,buffer[i]);
	}

	//convergence
	tMatrix aux(rows,columns);
	for(i=0;i<iterationsForConvergence;i++)
	{
		aux = 0;
		for(j=0;j<nCoefficients;j++)
			// aux = aux + coefficients[j]*buffer[(i+nCoefficientes-1-j) % nCoefficientes];
			Util::Add(aux,buffer[(i+nCoefficients-1-j) % nCoefficients],aux,1.0,coefficients[j]);

		noise = *(new tMatrix(randomGenerator->randnArray(rows*columns,noiseMean,noiseVariance),rows,columns));

		// buffer[i % nCoefficients] = aux + noise;
		Util::Add(aux,noise,buffer[i % nCoefficients]);

// 		cout << buffer[i % nCoefficients] << endl;
	}

	// the index of the next matrix is gonna be returned is kept
	iNextMatrix = i;
}

ARprocess::~ARprocess()
{
}

tMatrix ARprocess::NextMatrix()
{
// 	Matrix resultado = MatrixBuilder.CreateMatrix(nFilas,nColumnas);
//
// 	resultado.Clear();
// 	for(int j=0;j<nCoeficientes;j++)
// 		resultado = resultado + buffer[(iSiguienteMatriz+nCoeficientes-1-j) % (nCoeficientes)];
// 	resultado = resultado + RandomUtil.MatrizGaussiana(nFilas,nColumnas,mediaRuido,varianzaRuido);
// 	buffer[iSiguienteMatriz%nCoeficientes] = resultado;
// 	return resultado;
	tMatrix aux(rows,columns);

	aux = 0;
	for(int j=0;j<nCoefficients;j++)
		// aux = aux + coefficients[j]*buffer[(i+nCoefficientes-1-j) % nCoefficientes];
		Util::Add(aux,buffer[(iNextMatrix+nCoefficients-1-j) % nCoefficients],aux,1.0,coefficients[j]);

	tMatrix noise(randomGenerator->randnArray(rows*columns,noiseMean,noiseVariance),rows,columns);

	// buffer[i % nCoefficients] = aux + noise;
	Util::Add(aux,noise,buffer[iNextMatrix % nCoefficients]);

	return buffer[iNextMatrix % nCoefficients];
}


