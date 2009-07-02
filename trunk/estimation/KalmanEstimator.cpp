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

// #define DEBUG

// KalmanEstimator::KalmanEstimator(const tMatrix &initialEstimation,int N,double ARcoefficient,double ARvariance)
//  : ChannelMatrixEstimator(initialEstimation,N),_nChannelCoefficients(_L*_Nm)
// {
// 	tMatrix R = LaGenMatDouble::eye(_nChannelCoefficients);
// 	R *= ARcoefficient;
// 	tMatrix stateEquationCovariance = LaGenMatDouble::eye(_nChannelCoefficients);
// 	stateEquationCovariance *= ARvariance;
// 	tVector initialMeanVector = Util::ToVector(initialEstimation,rowwise);
// 	tMatrix initialCovariance = LaGenMatDouble::eye(_nChannelCoefficients);
// 
// 	_kalmanFilter = new KalmanFilter(R,stateEquationCovariance,initialMeanVector,initialCovariance);
// }

KalmanEstimator::KalmanEstimator(const tMatrix &initialEstimation,const tMatrix &variances,int N,vector<double> ARcoefficients,double ARvariance): ChannelMatrixEstimator(initialEstimation,N),_nChannelCoefficients(_L*_Nm)
{
	tMatrix R = LaGenMatDouble::eye(_nChannelCoefficients);
	R *= ARcoefficients[0];
	tMatrix stateEquationCovariance = LaGenMatDouble::eye(_nChannelCoefficients);
	stateEquationCovariance *= ARvariance;
	tVector initialMeanVector = Util::ToVector(initialEstimation,rowwise);
    tMatrix initialCovariance = LaGenMatDouble::from_diag(Util::ToVector(variances,rowwise));

	_kalmanFilter = new KalmanFilter(R,stateEquationCovariance,initialMeanVector,initialCovariance);
}

KalmanEstimator::KalmanEstimator(const KalmanEstimator &kalmanEstimator):ChannelMatrixEstimator(kalmanEstimator),_kalmanFilter(new KalmanFilter(*(kalmanEstimator._kalmanFilter))),_nChannelCoefficients(kalmanEstimator._nChannelCoefficients)
{
}

KalmanEstimator::~KalmanEstimator()
{
	delete _kalmanFilter;
}

tMatrix KalmanEstimator::nextMatrix(const tVector &observations,const tMatrix &symbolsMatrix,double noiseVariance)
{
	if(observations.size()!=_L || symbolsMatrix.rows()*symbolsMatrix.cols()!=_Nm)
		throw RuntimeException("KalmanEstimator::NextMatrix: observations vector length or symbols matrix length are wrong.");

	tMatrix observationEquationCovariance = LaGenMatDouble::eye(_L);
	observationEquationCovariance *= noiseVariance;
	_kalmanFilter->Step(BuildFfromSymbolsMatrix(Util::ToVector(symbolsMatrix,columnwise)),observations,observationEquationCovariance);

	_lastEstimatedChannelMatrix = Util::ToMatrix(_kalmanFilter->FilteredMean(),rowwise,_L);

	return  _lastEstimatedChannelMatrix;
}

tMatrix KalmanEstimator::BuildFfromSymbolsMatrix(const tVector &symbolsVector)
{
	int i,j;
	tMatrix res = LaGenMatDouble::zeros(_L,_L*_Nm);

	if(symbolsVector.size()!=_Nm)
		throw RuntimeException("KalmanEstimator::BuildFfromSymbolsMatrix: the number of elements of the received symbols vector is wrong.");

	// stacks the symbols inside symbolsMatrix to construct F
	for(i=0;i<_L;i++)
		for(j=0;j<_Nm;j++)
			res(i,i*_Nm+j) = symbolsVector(j);
	return res;
}

double KalmanEstimator::likelihood(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance)
{
	if(observations.size()!=_L || (symbolsMatrix.rows()*symbolsMatrix.cols())!=_Nm)
		throw RuntimeException("KalmanEstimator::likelihood: observations vector length or symbols matrix length are wrong.");

	// pivots vector needed for factorizations
	tLongIntVector piv(_nChannelCoefficients);

	tMatrix invPredictiveCovariance = _kalmanFilter->PredictiveCovariance();
	LUFactorizeIP(invPredictiveCovariance,piv);
    
	// detPredictiveCovariance = det(_kalmanFilter->PredictiveCovariance())
	double detPredictiveCovariance = 1.0;
	for(int i=0;i<_nChannelCoefficients;i++)
		detPredictiveCovariance *= invPredictiveCovariance(i,i);

	// invPredictiveCovariance = inv(_kalmanFilter->PredictiveCovariance())
	LaLUInverseIP(invPredictiveCovariance,piv);

	tMatrix F = BuildFfromSymbolsMatrix(Util::ToVector(symbolsMatrix,columnwise));

	tMatrix B = invPredictiveCovariance;
	// B = invPredictiveCovariance + (1.0/noiseVariance) F' * F
	Blas_Mat_Trans_Mat_Mult(F,F,B,1.0/noiseVariance,1.0);


	// B = inv(B) <------------------------------------------------------------------------
	LUFactorizeIP(B,piv);
	LaLUInverseIP(B,piv);

	tVector invPredictiveCovariancePredictiveMean(_nChannelCoefficients);

	// invPredictiveCovariancePredictiveMean = invPredictiveCovariance * _kalmanFilter->PredictiveMean()
	Blas_Mat_Vec_Mult(invPredictiveCovariance,_kalmanFilter->PredictiveMean(),invPredictiveCovariancePredictiveMean);

	tVector auxAuxArgExp = invPredictiveCovariancePredictiveMean;
	// auxAuxArgExp = invPredictiveCovariancePredictiveMean + (1.0/noiseVariance) F' * observations
	Blas_Mat_Trans_Vec_Mult(F,observations,auxAuxArgExp,1.0/noiseVariance,1.0);

	tVector auxAuxArgExpInvB(_nChannelCoefficients);

	// auxAuxArgExpInvB = B' * auxAuxArgExp = auxAuxArgExp * B
	Blas_Mat_Trans_Vec_Mult(B,auxAuxArgExp,auxAuxArgExpInvB);

	// _auxArgExp = auxAuxArgExpInvB . auxAuxArgExp
	double auxArgExp = Blas_Dot_Prod(auxAuxArgExpInvB,auxAuxArgExp);

	tVector observationsNoiseCovariance = observations;
// 	observationsNoiseCovariance *= noiseVariance;
    observationsNoiseCovariance *= 1.0/noiseVariance; // <--------------------------------------------------------------------------------------

	// _observationsNoiseCovarianceObservations = observationsNoiseCovariance . observations
	double observationsNoiseCovarianceObservations = Blas_Dot_Prod(observationsNoiseCovariance,observations);

	// predictiveMeanInvPredictiveCovariancePredictiveMean = _kalmanFilter->PredictiveMean() . invPredictiveCovariancePredictiveMean
	double predictiveMeanInvPredictiveCovariancePredictiveMean = Blas_Dot_Prod(_kalmanFilter->PredictiveMean(),invPredictiveCovariancePredictiveMean);

	double argExp = -0.5*(observationsNoiseCovarianceObservations + predictiveMeanInvPredictiveCovariancePredictiveMean - auxArgExp);

	//detInvB = det(B) (recall B = inv(B))
	LUFactorizeIP(B,piv);
	double detInvB = 1.0;
	for(int i=0;i<_nChannelCoefficients;i++)
		detInvB *= B(i,i);

	return sqrt(fabs(detInvB))/(pow(2*M_PI*noiseVariance,_L/2)*sqrt(fabs(detPredictiveCovariance)))*exp(argExp);
}

KalmanEstimator *KalmanEstimator::Clone() const
{
	// it relies on copy constructor
	return new KalmanEstimator(*this);
}

tMatrix KalmanEstimator::SampleFromPredictive()
{
	tVector predictiveMean = _kalmanFilter->PredictiveMean();
	tMatrix predictiveCovariance = _kalmanFilter->PredictiveCovariance();

	return Util::ToMatrix(StatUtil::RandMatrix(predictiveMean,predictiveCovariance),rowwise,_L);
}

void KalmanEstimator::setFirstEstimatedChannelMatrix(const tMatrix &matrix)
{
    ChannelMatrixEstimator::setFirstEstimatedChannelMatrix(matrix);

    _kalmanFilter->SetFilteredMean(Util::ToVector(matrix,columnwise));
}
