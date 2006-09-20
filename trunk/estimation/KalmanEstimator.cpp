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

KalmanEstimator::KalmanEstimator(tMatrix initialEstimation,double ARcoefficient,double ARvariance)
 : ChannelMatrixEstimator(initialEstimation),_nChannelCoefficients(_L*_Nm),_identityL(LaGenMatDouble::eye(_L)),
//auxiliary variables initialization
_F(LaGenMatDouble::zeros(_L,_L*_Nm)),_piv(_nChannelCoefficients),_FtransInvNoiseCovariance(_nChannelCoefficients,_L),_B(_nChannelCoefficients,_nChannelCoefficients),_invPredictiveCovariancePredictiveMean(_nChannelCoefficients),_auxAuxArgExp(_nChannelCoefficients),_auxAuxArgExpInvB(_nChannelCoefficients),_observationsNoiseCovariance(_L)
{
	tMatrix R = LaGenMatDouble::eye(_nChannelCoefficients);
	R *= ARcoefficient;
	tMatrix stateEquationCovariance = LaGenMatDouble::eye(_nChannelCoefficients);
	stateEquationCovariance *= ARvariance;
	tVector initialMeanVector = Util::ToVector(initialEstimation,rowwise);
	tMatrix initialCovariance = LaGenMatDouble::eye(_nChannelCoefficients);

	_kalmanFilter = new KalmanFilter(R,stateEquationCovariance,initialMeanVector,initialCovariance,_L);
}

KalmanEstimator::KalmanEstimator(const KalmanEstimator &kalmanEstimator):ChannelMatrixEstimator(kalmanEstimator),_kalmanFilter(new KalmanFilter(*(kalmanEstimator._kalmanFilter))),_nChannelCoefficients(kalmanEstimator._nChannelCoefficients),_identityL(kalmanEstimator._identityL),_F(LaGenMatDouble::zeros(_L,_L*_Nm)),_piv(_nChannelCoefficients),_FtransInvNoiseCovariance(_nChannelCoefficients,_L),_B(_nChannelCoefficients,_nChannelCoefficients),_invPredictiveCovariancePredictiveMean(_nChannelCoefficients),_auxAuxArgExp(_nChannelCoefficients),_auxAuxArgExpInvB(_nChannelCoefficients),_observationsNoiseCovariance(_L)
{
}

KalmanEstimator::~KalmanEstimator()
{
	delete _kalmanFilter;
}

tMatrix KalmanEstimator::NextMatrix(const tVector &observations,const tMatrix &symbolsMatrix,double noiseVariance)
{
// 	cout << "obser es " << observations.size() << " L es " << _L << " symbolsMatrix.rows()*symbolsMatrix.cols() " << symbolsMatrix.rows()*symbolsMatrix.cols() << " _Nm " << _Nm << endl;
	if(observations.size()!=_L || symbolsMatrix.rows()*symbolsMatrix.cols()!=_Nm)
		throw RuntimeException("KalmanEstimator::NextMatrix: observations vector length or symbols matrix length are wrong.");

	FillFfromSymbolsMatrix(symbolsMatrix);
	tMatrix observationEquationCovariance = LaGenMatDouble::eye(_L);
	observationEquationCovariance *= noiseVariance;
	_kalmanFilter->Step(_F,observations,observationEquationCovariance);

	_lastEstimatedChannelMatrix = Util::ToMatrix(_kalmanFilter->FilteredMean(),rowwise,_L);
		
	return  _lastEstimatedChannelMatrix;
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

double KalmanEstimator::Likelihood(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance)
{
// 	cout << "obser es " << observations.size() << " L es " << _L << " symbolsMatrix.rows()*symbolsMatrix.cols() " << symbolsMatrix.rows()*symbolsMatrix.cols() << " _Nm " << _Nm << endl;
	if(observations.size()!=_L || (symbolsMatrix.rows()*symbolsMatrix.cols())!=_Nm)
		throw RuntimeException("KalmanEstimator::Likelihood: observations vector length or symbols matrix length are wrong.");

	tMatrix noiseCovariance = LaGenMatDouble::eye(_L);
	noiseCovariance *= noiseVariance;
	tMatrix invNoiseCovariance = LaGenMatDouble::eye(_L);
	invNoiseCovariance *= (1/noiseVariance);

	tMatrix invPredictiveCovariance = _kalmanFilter->PredictiveCovariance();
	LUFactorizeIP(invPredictiveCovariance,_piv);	
	// detPredictiveCovariance = det(_kalmanFilter->PredictiveCovariance())
	double detPredictiveCovariance = 1.0;
	for(int i=0;i<_nChannelCoefficients;i++)
		detPredictiveCovariance *= invPredictiveCovariance(i,i);

	// invPredictiveCovariance = inv(_kalmanFilter->PredictiveCovariance())
	LaLUInverseIP(invPredictiveCovariance,_piv);

	FillFfromSymbolsMatrix(symbolsMatrix);

	// _FtransInvNoiseCovariance = _F' * invNoiseCovariance
	Blas_Mat_Trans_Mat_Mult(_F,invNoiseCovariance,_FtransInvNoiseCovariance);

	// _B = _FtransInvNoiseCovariance * _F
	Blas_Mat_Mat_Mult(_FtransInvNoiseCovariance,_F,_B);

	// _B = _B + invPredictiveCovariance
	Util::Add(_B,invPredictiveCovariance,_B);

	// invB = inv(_B)
	tMatrix invB = _B;
	LUFactorizeIP(invB,_piv);	
	LaLUInverseIP(invB,_piv);

	// _invPredictiveCovariancePredictiveMean = invPredictiveCovariance * _kalmanFilter->PredictiveMean()
	Blas_Mat_Vec_Mult(invPredictiveCovariance,_kalmanFilter->PredictiveMean(),_invPredictiveCovariancePredictiveMean);

	// _auxAuxArgExp = _FtransInvNoiseCovariance * observations
	Blas_Mat_Vec_Mult(_FtransInvNoiseCovariance,observations,_auxAuxArgExp);

	// _auxAuxArgExp = _auxAuxArgExp + _invPredictiveCovariancePredictiveMean
	Util::Add(_auxAuxArgExp,_invPredictiveCovariancePredictiveMean,_auxAuxArgExp);

	// _auxAuxArgExpInvB = invB' * _auxAuxArgExp = _auxAuxArgExp * invB
	Blas_Mat_Trans_Vec_Mult(invB,_auxAuxArgExp,_auxAuxArgExpInvB);

	// _auxArgExp = _auxAuxArgExpInvB . _auxAuxArgExp
	double auxArgExp = Blas_Dot_Prod(_auxAuxArgExpInvB,_auxAuxArgExp);

// _observationsNoiseCovariance = noiseCovariance' * observations = observations * noiseCovariance
	Blas_Mat_Trans_Vec_Mult(noiseCovariance,observations,_observationsNoiseCovariance);

	// _observationsNoiseCovarianceObservations = _observationsNoiseCovariance . observations
	double observationsNoiseCovarianceObservations = Blas_Dot_Prod(_observationsNoiseCovariance,observations);

	// predictiveMeanInvPredictiveCovariancePredictiveMean = _kalmanFilter->PredictiveMean() . _invPredictiveCovariancePredictiveMean
	double predictiveMeanInvPredictiveCovariancePredictiveMean = Blas_Dot_Prod(_kalmanFilter->PredictiveMean(),_invPredictiveCovariancePredictiveMean);

	double argExp = -0.5*(observationsNoiseCovarianceObservations + predictiveMeanInvPredictiveCovariancePredictiveMean - auxArgExp);

	//detInvB = inv(invB)
	LUFactorizeIP(invB,_piv);
	double detInvB = 1.0;
	for(int i=0;i<_nChannelCoefficients;i++)
		detInvB *= invB(i,i);

	return sqrt(fabs(detInvB))/(pow(2*M_PI*noiseVariance,_L/2)*sqrt(fabs(detPredictiveCovariance)))*exp(argExp);
}

KalmanEstimator *KalmanEstimator::Clone()
{
	// it relies on copy constructor
	return new KalmanEstimator(*this);
}
