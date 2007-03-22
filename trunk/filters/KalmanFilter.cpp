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
#include "KalmanFilter.h"

KalmanFilter::KalmanFilter(const tMatrix &R,const tMatrix &stateEquationCovariance,const tVector &initialMean,const tMatrix &initialCovariance,int observationVectorLength):
_R(R),_stateEquationCovariance(stateEquationCovariance),_nElementsToEstimate(initialMean.size()),_predictiveMean(R.rows()),_filteredMean(initialMean),_predictiveCovariance(R.rows(),R.rows()),_filteredCovariance(initialCovariance)
{
	if(R.rows()!=_nElementsToEstimate || _nElementsToEstimate!=R.cols())
		throw RuntimeException("KalmanFilter::KalmanFilter: matrices R and F dimensions are not coherent with those of the vector to be estimated.");

	if(initialMean.size()!=initialCovariance.rows() || initialMean.size()!=initialCovariance.cols())
		throw RuntimeException("KalmanFilter::KalmanFilter: the number of rows and columns of the covariance must be the number of rows of the mean.");

	// _predictiveMean = _R*_filteredMean
	Blas_Mat_Vec_Mult(_R,_filteredMean,_predictiveMean);

	tMatrix RfilteredCovariance(_nElementsToEstimate,_nElementsToEstimate);
	// RfilteredCovariance = _R*_filteredCovariance
	Blas_Mat_Mat_Mult(_R,_filteredCovariance,RfilteredCovariance);

    _predictiveCovariance = _stateEquationCovariance;

    // _predictiveCovariance = _predictiveCovariance + RfilteredCovariance*_R'
    // note that _predictiveCovariance is initalized to _stateEquationCovariance
    Blas_Mat_Mat_Trans_Mult(RfilteredCovariance,_R,_predictiveCovariance,1.0,1.0);
}

void KalmanFilter::Step(const tMatrix &F,const tVector &observation,const tMatrix &observationEquationCovariance)
{
	if(F.cols()!=_nElementsToEstimate || F.rows()!=observation.size()/* || observation.size()!=_observationVectorLength*/)
		throw RuntimeException("The matrix F or observation vector dimensions are wrong.");

	tMatrix predictiveCovarianceFtrans(_nElementsToEstimate,observation.size());
	// predictiveCovarianceFtrans = _predictiveCovariance * F'
	Blas_Mat_Mat_Trans_Mult(_predictiveCovariance,F,predictiveCovarianceFtrans);

    tMatrix auxMatrix = observationEquationCovariance;
    // auxMatrix = auxMatrix + F * predictiveCovarianceFtrans
    Blas_Mat_Mat_Mult(F,predictiveCovarianceFtrans,auxMatrix,1.0,1.0);

	tLongIntVector piv(observation.size());
	// auxMatrix = inverse(auxMatrix)
	LUFactorizeIP(auxMatrix,piv);
	LaLUInverseIP(auxMatrix,piv);

	tMatrix KalmanGain(_nElementsToEstimate,observation.size());
	// KalmanGain = predictiveCovarianceFtrans * auxMatrix
	Blas_Mat_Mat_Mult(predictiveCovarianceFtrans,auxMatrix,KalmanGain);

    tVector auxVector = observation;
    // auxVector = auxVector - F * _predictiveMean
    Blas_Mat_Vec_Mult(F,_predictiveMean,auxVector,-1.0,1.0);

    _filteredMean = _predictiveMean;
    // _filteredMean = _filteredMean + KalmanGain * auxVector
    Blas_Mat_Vec_Mult(KalmanGain,auxVector,_filteredMean,1.0,1.0);

	tMatrix FpredictiveCovariance(observation.size(),_nElementsToEstimate);
	// FpredictiveCovariance = F * _predictiveCovariance
	Blas_Mat_Mat_Mult(F,_predictiveCovariance,FpredictiveCovariance);

    _filteredCovariance = _predictiveCovariance;
    // _filteredCovariance = _filteredCovariance - KalmanGain * FpredictiveCovariance
    Blas_Mat_Mat_Mult(KalmanGain,FpredictiveCovariance,_filteredCovariance,-1.0,1.0);

	// **************** PREDICTIVES **********************

	// _predictiveMean = _R * _filteredMean
	Blas_Mat_Vec_Mult(_R,_filteredMean,_predictiveMean);

	tMatrix RfilteredCovariance(_nElementsToEstimate,_nElementsToEstimate);
	// RfilteredCovariance = _R*_filteredCovariance
	Blas_Mat_Mat_Mult(_R,_filteredCovariance,RfilteredCovariance);

    _predictiveCovariance = _stateEquationCovariance;

    // _predictiveCovariance = _predictiveCovariance + RfilteredCovariance*_R'
    // note that _predictiveCovariance is initalized to _stateEquationCovariance
    Blas_Mat_Mat_Trans_Mult(RfilteredCovariance,_R,_predictiveCovariance,1.0,1.0);
}

void KalmanFilter::SetFilteredMean(const tVector &filteredMean)
{
	if(filteredMean.size()!=_filteredMean.size())
		throw RuntimeException("KalmanFilter::SetFilteredMean: the size of the received vector is wrong.");

	_filteredMean = filteredMean;

	// _predictiveMean = _R*_filteredMean
	Blas_Mat_Vec_Mult(_R,_filteredMean,_predictiveMean);
}
