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
// 	cout << "F antes" << _F << endl;
	for(i=0;i<_L;i++)
		for(j=0;j<_Nm;j++)
		{
			_F(i,i*_Nm+j) = symbolsMatrix(j%symbolsMatrix.rows(),j/symbolsMatrix.rows());
// 			cout << "modifico el elemento (" << i << "," << (i*_Nm+j) << ")" << endl;
		}
// 	cout << "F despues" << _F << endl;
}

double KalmanEstimator::Verosimilitud(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance)
{
	if(observations.size()!=_L || (symbolsMatrix.rows()*symbolsMatrix.cols())!=_Nm)
		throw RuntimeException("observations vector length or symbols matrix length are wrong.");

	tMatrix noiseCovariance = LaGenMatDouble::eye(_L);
	noiseCovariance *= noiseVariance;
	tMatrix invNoiseCovariance = LaGenMatDouble::eye(_L);
	invNoiseCovariance *= (1/noiseVariance);

	// _invPredictiveCovariance = inv(_kalmanFilter->PredictiveCovariance())
	tMatrix _invPredictiveCovariance = _kalmanFilter->PredictiveCovariance();
	tLongIntVector _piv(_nChannelCoefficients);
	LUFactorizeIP(_invPredictiveCovariance,_piv);
	LaLUInverseIP(_invPredictiveCovariance,_piv);

	FillFfromSymbolsMatrix(symbolsMatrix);

	tMatrix _FtransInvNoiseCovariance(_nChannelCoefficients,_L);
	// _FtransInvNoiseCovariance = _F' * invNoiseCovariance
	Blas_Mat_Trans_Mat_Mult(_F,invNoiseCovariance,_FtransInvNoiseCovariance);

	tMatrix _B(_nChannelCoefficients,_nChannelCoefficients);
	// _B = _FtransInvNoiseCovariance * _F
	Blas_Mat_Mat_Mult(_FtransInvNoiseCovariance,_F,_B);

	// _B = _B + _invPredictiveCovariance
	Util::Add(_B,_invPredictiveCovariance,_B);

	// invB = inv(_B)
	tMatrix invB = _B;
	LUFactorizeIP(invB,_piv);
	LaLUInverseIP(invB,_piv);

	tVector _invPredictiveCovariancePredictiveMean(_nChannelCoefficients);
	// _invPredictiveCovariancePredictiveMean = _invPredictiveCovariance * _kalmanFilter->PredictiveMean()
	Blas_Mat_Vec_Mult(_invPredictiveCovariance,_kalmanFilter->PredictiveMean(),_invPredictiveCovariancePredictiveMean);

	tVector _auxAuxArgExp(_nChannelCoefficients);
	// _auxAuxArgExp = _FtransInvNoiseCovariance * observations
	Blas_Mat_Vec_Mult(_FtransInvNoiseCovariance,observations,_auxAuxArgExp);

	// _auxAuxArgExp = _auxAuxArgExp + _invPredictiveCovariancePredictiveMean
	Util::Add(_auxAuxArgExp,_invPredictiveCovariancePredictiveMean,_auxAuxArgExp);

	tVector _auxAuxArgExpInvB(_nChannelCoefficients);
	// _auxAuxArgExpInvB = invB' * _auxAuxArgExp = _auxAuxArgExp * invB
	Blas_Mat_Trans_Vec_Mult(invB,_auxAuxArgExp,_auxAuxArgExpInvB);

}

