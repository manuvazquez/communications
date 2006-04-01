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
#include "KalmanEstimator.h"

KalmanEstimator::KalmanEstimator(double ARcoefficient,double ARvariance,tMatrix &initialMeanMatrix)
 : ChannelMatrixEstimator(initialMeanMatrix.rows(),initialMeanMatrix.cols()),_nChannelCoefficients(_L*_Nm),_identityL(LaGenMatDouble::eye(_L)),_F(_L,_L*_Nm)
{
	tMatrix R = LaGenMatDouble::eye(_nChannelCoefficients);
	R *= ARcoefficient;
	tMatrix stateEquationCovariance = LaGenMatDouble::eye(_nChannelCoefficients);
	stateEquationCovariance *= ARvariance;
	tVector initialMeanVector = Util::ToVector(initialMeanMatrix,rowwise);
	tMatrix initialCovariance = LaGenMatDouble::eye(_nChannelCoefficients);

	_kalmanFilter = new KalmanFilter(R,stateEquationCovariance,initialMeanVector,initialCovariance,_L);
	_F = 0.0;
// 	_identityL = LaGenMatDouble::eye(_L);
}

tMatrix KalmanEstimator::NextMatrix(const tVector &observations,const tMatrix &symbolsMatrix,double noiseVariance)
{
	if(observations.size()!=_L || symbolsMatrix.rows()*symbolsMatrix.cols()!=_Nm)
		throw RuntimeException("observations vector length or symbols matrix length are wrong.");

	FillFfromSymbolsMatrix(symbolsMatrix);
	tMatrix observationEquationCovariance = LaGenMatDouble::eye(_L);
	observationEquationCovariance *= noiseVariance;
	_kalmanFilter->Step(_F,observations,observationEquationCovariance);
	return Util::ToMatrix(_kalmanFilter->FilteredMean(),rowwise,_L);
}

void KalmanEstimator::FillFfromSymbolsMatrix(const tMatrix &symbolsMatrix)
{
	int i,j,nSymbols;

	nSymbols = symbolsMatrix.rows()*symbolsMatrix.cols();
	if(nSymbols!=_Nm)
		throw RuntimeException("The number of elements of the received symbols matrix is wrong.");

	// stacks the symbols inside symbolsMatrix to construct F
	for(i=0;i<_L;i++)
		for(j=0;j<_nChannelCoefficients;j++)
			_F(i,i*_nChannelCoefficients+j) = symbolsMatrix(j/symbolsMatrix.cols(),j%symbolsMatrix.cols());
}

