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
_R(R),_stateEquationCovariance(stateEquationCovariance),_nElementsToEstimate(initialMean.size()),_predictiveMean(R.rows()),_filteredMean(initialMean),_predictiveCovariance(R.rows(),R.rows()),_filteredCovariance(initialCovariance),_observationVectorLength(observationVectorLength),
// auxiliary variables
_predictiveCovarianceFtrans(_nElementsToEstimate,_observationVectorLength)/*,_auxMatrix(_observationVectorLength,_observationVectorLength)*/,_KalmanGain(_nElementsToEstimate,_observationVectorLength),_FpredictiveCovariance(_observationVectorLength,_nElementsToEstimate)/*,_KalmanGainFpredictiveCovariance(_nElementsToEstimate,_nElementsToEstimate)*//*,_predictiveCovarianceAux(_nElementsToEstimate,_nElementsToEstimate)*/,_RfilteredCovariance(_nElementsToEstimate,_nElementsToEstimate)/*,_RfilteredCovarianceRtrans(_nElementsToEstimate,_nElementsToEstimate)*//*,_auxVector(_observationVectorLength)*//*,_KalmanGainByNotPredicted(_nElementsToEstimate)*/,_piv(_observationVectorLength)
{
	if(R.rows()!=_nElementsToEstimate || _nElementsToEstimate!=R.cols())
		throw RuntimeException("KalmanFilter::KalmanFilter: matrices R and F dimensions are not coherent with those of the vector to be estimated.");

	if(initialMean.size()!=initialCovariance.rows() || initialMean.size()!=initialCovariance.cols())
		throw RuntimeException("KalmanFilter::KalmanFilter: the number of rows and columns of the covariance must be the number of rows of the mean.");

	// _predictiveMean = _R*_filteredMean
	Blas_Mat_Vec_Mult(_R,_filteredMean,_predictiveMean);

	// _RfilteredCovariance = _R*_filteredCovariance
	Blas_Mat_Mat_Mult(_R,_filteredCovariance,_RfilteredCovariance);

// 	// _RfilteredCovarianceRtrans = _RfilteredCovariance*_R'
// 	Blas_Mat_Mat_Trans_Mult(_RfilteredCovariance,_R,_RfilteredCovarianceRtrans);
//
// 	// _predictiveCovariance = _RfilteredCovarianceRtrans + _stateEquationCovariance;
// 	Util::Add(_RfilteredCovarianceRtrans,_stateEquationCovariance,_predictiveCovariance);

    _predictiveCovariance = _stateEquationCovariance;

    // _predictiveCovariance = _predictiveCovariance + _RfilteredCovariance*_R'
    // note that _predictiveCovariance is initalized to _stateEquationCovariance
    Blas_Mat_Mat_Trans_Mult(_RfilteredCovariance,_R,_predictiveCovariance,1.0,1.0);

}

void KalmanFilter::Step(const tMatrix &F,const tVector &observation,const tMatrix &observationEquationCovariance)
{
	if(F.cols()!=_nElementsToEstimate || F.rows()!=observation.size() || observation.size()!=_observationVectorLength)
		throw RuntimeException("The matrix F or observation vector dimensions are wrong.");

	// _predictiveCovarianceFtrans = _predictiveCovariance * F'
	Blas_Mat_Mat_Trans_Mult(_predictiveCovariance,F,_predictiveCovarianceFtrans);

// 	// _auxMatrix = F * _predictiveCovarianceFtrans
// 	Blas_Mat_Mat_Mult(F,_predictiveCovarianceFtrans,_auxMatrix);
//
// 	// _auxMatrix = _auxMatrix + observationEquationCovariance
// 	Util::Add(_auxMatrix,observationEquationCovariance,_auxMatrix);


    tMatrix auxMatrix = observationEquationCovariance;
    // auxMatrix = auxMatrix + F * _predictiveCovarianceFtrans
    Blas_Mat_Mat_Mult(F,_predictiveCovarianceFtrans,auxMatrix,1.0,1.0);




	// auxMatrix = inverse(auxMatrix)
	LUFactorizeIP(auxMatrix,_piv);
	LaLUInverseIP(auxMatrix,_piv);

	// _KalmanGain = _predictiveCovarianceFtrans * auxMatrix
	Blas_Mat_Mat_Mult(_predictiveCovarianceFtrans,auxMatrix,_KalmanGain);

// 	// _auxVector = F * _predictiveMean
// 	Blas_Mat_Vec_Mult(F,_predictiveMean,_auxVector);
//
// 	// _auxVector = observation - _auxVector
// 	Util::Add(observation,_auxVector,_auxVector,1.0,-1.0);


    tVector auxVector = observation;
    // auxVector = auxVector - F * _predictiveMean
    Blas_Mat_Vec_Mult(F,_predictiveMean,auxVector,-1.0,1.0);


// 	// 	_KalmanGainByNotPredicted = _KalmanGain * auxVector
// 	Blas_Mat_Vec_Mult(_KalmanGain,auxVector,_KalmanGainByNotPredicted);
//
// 	// _filteredMean = _predictiveMean + _KalmanGainByNotPredicted
// 	Util::Add(_predictiveMean,_KalmanGainByNotPredicted,_filteredMean);


    _filteredMean = _predictiveMean;
    // _filteredMean = _filteredMean + _KalmanGain * auxVector
    Blas_Mat_Vec_Mult(_KalmanGain,auxVector,_filteredMean,1.0,1.0);


	// _FpredictiveCovariance = F * _predictiveCovariance
	Blas_Mat_Mat_Mult(F,_predictiveCovariance,_FpredictiveCovariance);

// 	// _KalmanGainFpredictiveCovariance = _KalmanGain * _FpredictiveCovariance
// 	Blas_Mat_Mat_Mult(_KalmanGain,_FpredictiveCovariance,_KalmanGainFpredictiveCovariance);
//
// 	// _filteredCovariance = _predictiveCovariance - _KalmanGainFpredictiveCovariance
// 	Util::Add(_predictiveCovariance,_KalmanGainFpredictiveCovariance,_filteredCovariance,1.0,-1.0);


    _filteredCovariance = _predictiveCovariance;
    // _filteredCovariance = _filteredCovariance - _KalmanGain * _FpredictiveCovariance
    Blas_Mat_Mat_Mult(_KalmanGain,_FpredictiveCovariance,_filteredCovariance,-1.0,1.0);

	// **************** PREDICTIVES **********************

	// _predictiveMean = _R * _filteredMean
	Blas_Mat_Vec_Mult(_R,_filteredMean,_predictiveMean);

	// _RfilteredCovariance = _R*_filteredCovariance
	Blas_Mat_Mat_Mult(_R,_filteredCovariance,_RfilteredCovariance);

// 	// _RfilteredCovarianceRtrans = _R*_filteredCovariance*_R'
// 	Blas_Mat_Mat_Trans_Mult(_RfilteredCovariance,_R,_RfilteredCovarianceRtrans);
//
// 	// _predictiveCovariance = _RfilteredCovarianceRtrans + _stateEquationCovariance;
// 	Util::Add(_RfilteredCovarianceRtrans,_stateEquationCovariance,_predictiveCovariance);

    _predictiveCovariance = _stateEquationCovariance;

    // _predictiveCovariance = _predictiveCovariance + _RfilteredCovariance*_R'
    // note that _predictiveCovariance is initalized to _stateEquationCovariance
    Blas_Mat_Mat_Trans_Mult(_RfilteredCovariance,_R,_predictiveCovariance,1.0,1.0);
}
