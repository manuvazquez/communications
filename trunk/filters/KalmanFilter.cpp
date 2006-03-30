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

KalmanFilter::KalmanFilter(tMatrix R,tMatrix stateEquationCovariance,tVector initialMean, tMatrix initialCovariance,int observationVectorLength):
_nElementsToEstimate(initialMean.size()),_R(R),_stateEquationCovariance(stateEquationCovariance),_filteredMean(initialMean),_filteredCovariance(initialCovariance),_predictiveMean(R.rows()),_predictiveCovariance(R.rows(),R.rows()),_observationVectorLength(observationVectorLength)
{
	if(R.rows()!=_nElementsToEstimate || _nElementsToEstimate!=R.cols())
		throw RuntimeException("Matrices R and F dimensions are not coherent with those of the vector to be estimated.");
	
	// _predictiveMean = _R*_filteredMean
	Blas_Mat_Vec_Mult(_R,_filteredMean,_predictiveMean);

	tMatrix aux(_R.rows(),_R.rows());
	// aux = _R*_filteredCovariance
	Blas_Mat_Mat_Mult(_R,_filteredCovariance,aux);

	// aux = _R*_filteredCovariance*_R^T
	Blas_Mat_Mat_Trans_Mult(aux,_R,aux);

	// _predictiveCovariance = aux + _stateEquationCovariance;
	Util::Add(aux,_stateEquationCovariance,_predictiveCovariance);
	
	// memory for several auxiliar matrices is allocated
	_predictiveCovarianceF(_nElementsToEstimate,_observationVectorLength);
	_auxMatrix(_observationVectorLength,_observationVectorLength);
}


KalmanFilter::~KalmanFilter()
{
}

void KalmanFilter::Step(tMatrix F,tVector observation, tMatrix observationEquationCovariance)
{
	if(F.cols()!=_nElementsToEstimate || F.rows()!=observation.size() || observation.size()!=_observationVectorLength)
		throw RuntimeException("The matrix F or observation vector dimensions are wrong.");

// 	tVector auxVector(observation.size());
		
	// _predictiveCovarianceF = _predictiveCovariance * F'
	Blas_Mat_Mat_Trans_Mult(_predictiveCovariance,F,_predictiveCovarianceF);

	// _auxMatrix = F * _predictiveCovarianceF
	Blas_Mat_Mat_Mult(F,_predictiveCovarianceF,_auxMatrix);

	// _auxMatrix = _auxMatrix + observationEquationCovariance
	Util::Add(_auxMatrix,observationEquationCovariance,_auxMatrix);
	
// 	tLongIntVector pivotes(matrizAleatoria.size(0));
// 	LUFactorizeIP(matrizAleatoria,pivotes);
// 	LaLUInverseIP(matrizAleatoria,pivotes);
}
