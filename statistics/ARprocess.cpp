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

// #define DEBUG5

// #include <Eigen/Dense>
// #include <Eigen/LU>

// using namespace Eigen;

ARprocess::ARprocess(MatrixXd seed,vector<double> coefficients,double noiseVariance):_coefficients(coefficients),_noiseVariance(noiseVariance),_nCoefficients(coefficients.size()),_rows(seed.rows()),_columns(seed.cols())
{
    CommonConstructorsCode(seed);
}

ARprocess::ARprocess(MatrixXd seed,int order,double velocity,double carrierFrequency,double T):_rows(seed.rows()),_columns(seed.cols())
{
    _coefficients = parametersFromYuleWalker(order,velocity,carrierFrequency,T,_noiseVariance);

    _nCoefficients = _coefficients.size();

    CommonConstructorsCode(seed);
}

void ARprocess::CommonConstructorsCode(const MatrixXd &seed)
{
    int i,j;

    _iterationsForConvergence = 200;
    _noiseMean = 0.0;

    _buffer.resize(_nCoefficients);
    _buffer[0] = seed;

    //the buffer is filled
    for(i=1;i<_nCoefficients;i++)
        _buffer[i] = _buffer[i-1] + StatUtil::randnMatrix(_rows,_columns,_noiseMean,_noiseVariance);

    //convergence
    for(i=0;i<_iterationsForConvergence;i++)
    {
        MatrixXd aux = MatrixXd::Zero(_rows,_columns);
        for(j=0;j<_nCoefficients;j++)
            // aux = aux + _coefficients[j]*_buffer[(i+nCoefficientes-1-j) % nCoefficientes];
            aux += _coefficients[j]*_buffer[(i+_nCoefficients-1-j) % _nCoefficients];

        // _buffer[i % _nCoefficients] = aux + noise;
        _buffer[i % _nCoefficients] = aux + StatUtil::randnMatrix(_rows,_columns,_noiseMean,_noiseVariance);
    }

    // the index of the next matrix is gonna be returned is kept
    _iNextMatrix = i;
}

MatrixXd ARprocess::nextMatrix()
{
    MatrixXd aux = MatrixXd::Zero(_rows,_columns);
    for(int j=0;j<_nCoefficients;j++)
        // aux = aux + _coefficients[j]*_buffer[(i+nCoefficientes-1-j) % nCoefficientes];
        aux += _coefficients[j]*_buffer[(_iNextMatrix+_nCoefficients-1-j) % _nCoefficients];
    
    // _buffer[i % _nCoefficients] = aux + noise;
    _buffer[_iNextMatrix % _nCoefficients] = aux + StatUtil::randnMatrix(_rows,_columns,_noiseMean,_noiseVariance);

    return _buffer[_iNextMatrix++ % _nCoefficients];
}

vector<double> ARprocess::parametersFromYuleWalker(int order,double velocity,double carrierFrequency,double T,double &noiseVariance)
{
    const double c = 3e8;

	// (c/carrierFrequency) is the wavelength
    double dopplerFrequency = velocity/(c/carrierFrequency);

    double normDopplerFrequency = T*dopplerFrequency;

    MatrixXd autocorrelationsMatrix(order,order);
    VectorXd autocorrelationsVector(order);

	// for efficiency's sake, all needed correlations are computed here
    vector<double> autocorrelations(order+1);
    for(int i=0;i<=order;i++)
        autocorrelations[i] = j0(2.0*M_PI*normDopplerFrequency*double(i));

    for(int m=0;m<order;m++)
    {
        for(int k=0;k<order;k++)
            autocorrelationsMatrix(m,k) = autocorrelations[abs(m-k)];

        autocorrelationsVector(m) = autocorrelations[m+1];
    }

//     VectorXd coefficients;
//     autocorrelationsMatrix.lu().solve(autocorrelationsVector,&coefficients);
	
	PartialPivLU<MatrixXd> luDecomp(autocorrelationsMatrix);
	VectorXd coefficients = luDecomp.solve(autocorrelationsVector);

    // the variance of the noise will be computed in the loop...
    noiseVariance = autocorrelations[0];

    // ...besides, the eigen vector will be copied to a standard c++ vector
    vector<double> coefficientsCppVector(order);
    for(int i=0;i<order;i++)
    {
        noiseVariance -= coefficients(i)*autocorrelations[i+1];
        coefficientsCppVector[i] = coefficients(i);
    }

    // it takes into accout rounding errors
    if(noiseVariance<0.0)
        noiseVariance = 0.0;

    return coefficientsCppVector;
}
