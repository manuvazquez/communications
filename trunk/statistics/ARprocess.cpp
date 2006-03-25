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

ARprocess::ARprocess()
{
}

ARprocess::ARprocess(tMatrix seed,vector<double> coefficients,double noiseVariance)
{
	this->coefficients = coefficients;
	this->noiseVariance = noiseVariance;
	this->rows = seed.rows();
	this->columns = seed.cols();
	this->nCoefficients = coefficients.size();
	this->buffer = new tMatrix[nCoefficients];

// 			for(i=1;i<nCoeficientes;i++)
// 			{
// 				for(j=0;j<i;j++)
// 					buffer[i] = buffer[i] + buffer[j]*coeficientesAR[j];
// 				buffer[i] = buffer[i] + RandomUtil.MatrizGaussiana(nFilas,nColumnas,mediaRuido,varianzaRuido);
// 			}
// 			// ... y comienza el proceso iterativo
// 			Matrix aux = MatrixBuilder.CreateMatrix(nFilas,nColumnas);
// 			for(i=0;i<nIteracionesHastaConvergencia;i++)
// 			{
// 				aux.Clear();
// 				for(j=0;j<nCoeficientes;j++)
// 					aux = aux + buffer[(i+nCoeficientes-1-j) % (nCoeficientes)];
// 				buffer[i%nCoeficientes] = aux + RandomUtil.MatrizGaussiana(nFilas,nColumnas,mediaRuido,varianzaRuido);
// 				//System.Console.WriteLine("{0}\n",buffer[i%nCoeficientes]);
// 			}
//
// 			// se guarda el indice de ultima matriz generada
// 			iSiguienteMatriz = i;
	//the buffer is filled
	int i,j;
	buffer[0] = seed;
	for(i=1;i<nCoefficients;i++)
	{
		buffer[i] = *(new tMatrix(rows,columns));
		buffer[i] = 0;
		for(j=0;j<i;j++)
			Util::Add(buffer[i],buffer[j],buffer[i],1.0,coefficients[j]);
	}
}

ARprocess::~ARprocess()
{
}

// 	public class ProcesoAR
// 	{
// 	vector<double> coefficients;
// 	double noiseVariance;
// 	double noiseMean;
// 	int nCoefficients, rows, columns, nextMatrix;
// 	int iterationsForConvergence;
// 	tMatrix *buffer;
// //
// 		public ProcesoAR(Matrix matrizInicial,double[] coeficientesAR,double varianzaRuido)
// 		{
// 			this.coeficientesAR = coeficientesAR;
// 			this.varianzaRuido = varianzaRuido;
// 			nFilas = matrizInicial.Rows;
// 			nColumnas = matrizInicial.Columns;
//
// 			nCoeficientes = coeficientesAR.Length;
//
// 			// se incializa el buffer...
// 			buffer = new Matrix[nCoeficientes];
//
// 			// ... se llena...
// 			int i,j;
// 			buffer[0] = matrizInicial;
// 			for(i=1;i<nCoeficientes;i++)
// 			{
// 				for(j=0;j<i;j++)
// 					buffer[i] = buffer[i] + buffer[j]*coeficientesAR[j];
// 				buffer[i] = buffer[i] + RandomUtil.MatrizGaussiana(nFilas,nColumnas,mediaRuido,varianzaRuido);
// 			}
//
// 			// ... y comienza el proceso iterativo
// 			Matrix aux = MatrixBuilder.CreateMatrix(nFilas,nColumnas);
// 			for(i=0;i<nIteracionesHastaConvergencia;i++)
// 			{
// 				aux.Clear();
// 				for(j=0;j<nCoeficientes;j++)
// 					aux = aux + buffer[(i+nCoeficientes-1-j) % (nCoeficientes)];
// 				buffer[i%nCoeficientes] = aux + RandomUtil.MatrizGaussiana(nFilas,nColumnas,mediaRuido,varianzaRuido);
// 				//System.Console.WriteLine("{0}\n",buffer[i%nCoeficientes]);
// 			}
//
// 			// se guarda el indice de ultima matriz generada
// 			iSiguienteMatriz = i;
// 		}
