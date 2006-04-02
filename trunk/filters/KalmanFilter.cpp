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
_nElementsToEstimate(initialMean.size()),_R(R),_stateEquationCovariance(stateEquationCovariance),_filteredMean(initialMean),_filteredCovariance(initialCovariance),_predictiveMean(R.rows()),_predictiveCovariance(R.rows(),R.rows()),_observationVectorLength(observationVectorLength),
// auxiliary variables
_RfilteredCovariance(_nElementsToEstimate,_nElementsToEstimate),_RfilteredCovarianceRtrans(_nElementsToEstimate,_nElementsToEstimate),_predictiveCovarianceFtrans(_nElementsToEstimate,_observationVectorLength),_auxMatrix(_observationVectorLength,_observationVectorLength),_KalmanGain(_nElementsToEstimate,_observationVectorLength),_auxVector(_observationVectorLength),_KalmanGainByNotPredicted(_nElementsToEstimate),_FpredictiveCovariance(_observationVectorLength,_nElementsToEstimate),_KalmanGainFpredictiveCovariance(_nElementsToEstimate,_nElementsToEstimate),_predictiveCovarianceAux(_nElementsToEstimate,_nElementsToEstimate),_piv(_observationVectorLength)
{
	if(R.rows()!=_nElementsToEstimate || _nElementsToEstimate!=R.cols())
		throw RuntimeException("Matrices R and F dimensions are not coherent with those of the vector to be estimated.");

	// _predictiveMean = _R*_filteredMean
	Blas_Mat_Vec_Mult(_R,_filteredMean,_predictiveMean);

	// _RfilteredCovariance = _R*_filteredCovariance
	Blas_Mat_Mat_Mult(_R,_filteredCovariance,_RfilteredCovariance);

	// _RfilteredCovarianceRtrans = _RfilteredCovariance*_R'
	Blas_Mat_Mat_Trans_Mult(_RfilteredCovariance,_R,_RfilteredCovarianceRtrans);

	// _predictiveCovariance = _RfilteredCovarianceRtrans + _stateEquationCovariance;
	Util::Add(_RfilteredCovarianceRtrans,_stateEquationCovariance,_predictiveCovariance);

// 	cout << "Media filtrada" << endl << _filteredMean << endl;
// 	cout << "Media predictiva" << endl << _predictiveMean << endl;
// 	cout << "Covarianza estado" << endl << _stateEquationCovariance << endl;
// 	cout << "R" << endl << _R << endl;
//
// 	char c;
// 	cin >> c;

// 	// memory for several auxiliar matrices is allocated
// 	_predictiveCovarianceFtrans(_nElementsToEstimate,_observationVectorLength);
// 	_auxMatrix(_observationVectorLength,_observationVectorLength);
// 	_KalmanGain(_nElementsToEstimate,_observationVectorLength);
// 	_auxVector(_observationVectorLength);
// 	_KalmanGainByNotPredicted(_nElementsToEstimate);
// 	_FpredictiveCovariance(_observationVectorLength,_nElementsToEstimate);
// 	_KalmanGainFpredictiveCovariance(_nElementsToEstimate,_nElementsToEstimate);
// 	_predictiveCovarianceAux(_nElementsToEstimate,_nElementsToEstimate);
}

void KalmanFilter::Step(tMatrix F,tVector observation, tMatrix observationEquationCovariance)
{
	if(F.cols()!=_nElementsToEstimate || F.rows()!=observation.size() || observation.size()!=_observationVectorLength)
		throw RuntimeException("The matrix F or observation vector dimensions are wrong.");

	// _predictiveCovarianceFtrans = _predictiveCovariance * F'
	Blas_Mat_Mat_Trans_Mult(_predictiveCovariance,F,_predictiveCovarianceFtrans);

	// _auxMatrix = F * _predictiveCovarianceFtrans
	Blas_Mat_Mat_Mult(F,_predictiveCovarianceFtrans,_auxMatrix);

	// _auxMatrix = _auxMatrix + observationEquationCovariance
	Util::Add(_auxMatrix,observationEquationCovariance,_auxMatrix);

	// _auxMatrix = inverse(_auxMatrix)
	LUFactorizeIP(_auxMatrix,_piv);
	LaLUInverseIP(_auxMatrix,_piv);

	// _KalmanGain = _predictiveCovarianceFtrans * _auxMatrix
	Blas_Mat_Mat_Mult(_predictiveCovarianceFtrans,_auxMatrix,_KalmanGain);

	// _auxVector = F * _predictiveMean
	Blas_Mat_Vec_Mult(F,_predictiveMean,_auxVector);

	// _auxVector = observation - _auxVector
	Util::Add(observation,_auxVector,_auxVector,1.0,-1.0);

	// 	_KalmanGainByNotPredicted = _KalmanGain * _auxVector
	Blas_Mat_Vec_Mult(_KalmanGain,_auxVector,_KalmanGainByNotPredicted);

	// _filteredMean = _predictiveMean * _KalmanGainByNotPredicted
	Util::Add(_predictiveMean,_KalmanGainByNotPredicted,_filteredMean);

	// _FpredictiveCovariance = F * _predictiveCovariance
	Blas_Mat_Mat_Mult(F,_predictiveCovariance,_FpredictiveCovariance);

	// _KalmanGainFpredictiveCovariance = _KalmanGain * _FpredictiveCovariance
	Blas_Mat_Mat_Mult(_KalmanGain,_FpredictiveCovariance,_KalmanGainFpredictiveCovariance);

	// _filteredCovariance = _predictiveCovariance - _KalmanGainFpredictiveCovariance
	Util::Add(_predictiveCovariance,_KalmanGainFpredictiveCovariance,_filteredCovariance,1.0,-1.0);

	// **************** PREDICTIVES **********************

	// _predictiveMean = _R * _filteredMean
	Blas_Mat_Vec_Mult(_R,_filteredMean,_predictiveMean);

	// _RfilteredCovariance = _R*_filteredCovariance
	Blas_Mat_Mat_Mult(_R,_filteredCovariance,_RfilteredCovariance);

	// _RfilteredCovarianceRtrans = _R*_filteredCovariance*_R'
	Blas_Mat_Mat_Trans_Mult(_RfilteredCovariance,_R,_RfilteredCovarianceRtrans);

	// _predictiveCovariance = _RfilteredCovarianceRtrans + _stateEquationCovariance;
	Util::Add(_RfilteredCovarianceRtrans,_stateEquationCovariance,_predictiveCovariance);
}

//   		public virtual FiltroKalman Clone()
//   		{
//   			FiltroKalman clon = this.MemberwiseClone() as FiltroKalman;
//
//   			//clon.R = R.Clone();
//   			//clon.covarEcuacionEstado = covarEcuacionEstado.Clone();
//   			clon.mediaPredictiva = mediaPredictiva.Clone();
//   			clon.mediaFiltrada = mediaFiltrada.Clone();
//   			clon.covarPredictiva = covarPredictiva.Clone();
//   			clon.covarFiltrada = covarFiltrada.Clone();
// //  			Console.WriteLine("Holita\n");
// //  			covarPredictiva[0,1] = -100;
// //  			Console.WriteLine("No clon\n{0}\nClon\n{1}\n",covarPredictiva,clon.CovarianzaPredictiva);
//   			return clon;
//   		}
