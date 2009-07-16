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

// #define PRINT_INFO

KalmanFilter::KalmanFilter(const tMatrix &R,const tMatrix &stateEquationCovariance,const tVector &initialMean,const tMatrix &initialCovariance):
_R(R),_stateEquationCovariance(stateEquationCovariance),_nElementsToEstimate(initialMean.size()),_predictiveMean(R.rows()),_predictiveCovariance(R.rows(),R.rows()){
    if(R.rows()!=_nElementsToEstimate || _nElementsToEstimate!=R.cols())
        throw RuntimeException("KalmanFilter::KalmanFilter: matrices R and F dimensions are not coherent with those of the vector to be estimated.");

    if(initialMean.size()!=initialCovariance.rows() || initialMean.size()!=initialCovariance.cols())
        throw RuntimeException("KalmanFilter::KalmanFilter: the number of rows and columns of the covariance must be the number of rows of the mean.");
    
    setFilteredMean(initialMean);
    setFilteredCovariance(initialCovariance);
}

void KalmanFilter::step(const tMatrix &F,const tVector &observation,const tMatrix &observationEquationCovariance)
{
    if(F.cols()!=_nElementsToEstimate || F.rows()!=observation.size())
        throw RuntimeException("The matrix F or observation vector dimensions are wrong.");

    tMatrix predictiveCovarianceFtrans(_nElementsToEstimate,observation.size());
    // predictiveCovarianceFtrans = _predictiveCovariance * F'
    Blas_Mat_Mat_Trans_Mult(_predictiveCovariance,F,predictiveCovarianceFtrans);

    tMatrix residualCovariance = observationEquationCovariance;
    // residualCovariance = residualCovariance + F * predictiveCovarianceFtrans
    Blas_Mat_Mat_Mult(F,predictiveCovarianceFtrans,residualCovariance,1.0,1.0);

    tLongIntVector piv(observation.size());
    // residualCovariance = inverse(residualCovariance)
    LUFactorizeIP(residualCovariance,piv);
    LaLUInverseIP(residualCovariance,piv);

    tMatrix KalmanGain(_nElementsToEstimate,observation.size());
    // KalmanGain = predictiveCovarianceFtrans * residualCovariance
    Blas_Mat_Mat_Mult(predictiveCovarianceFtrans,residualCovariance,KalmanGain);

    tVector measurementResidual = observation;
    // measurementResidual = measurementResidual - F * _predictiveMean
    Blas_Mat_Vec_Mult(F,_predictiveMean,measurementResidual,-1.0,1.0);
    
#ifdef PRINT_INFO
    cout << "KalmanFilter::step: F is:" << endl << F;
    cout << "KalmanFilter::step: predictiveCovarianceFtrans is:" << endl << predictiveCovarianceFtrans;
    cout << "KalmanFilter::step: predictive covariance is:" << endl << _predictiveCovariance;
    cout << "KalmanFilter::step: residual covariance is:" << endl << residualCovariance;
    cout << "KalmanFilter::step: Kalman gain is:" << endl << KalmanGain;
    cout << "KalmanFilter::step: innovation is:" << endl << measurementResidual;
    cout << "KalmanFilter::step: filtered mean before update:" << endl << _filteredMean;
#endif    

    _filteredMean = _predictiveMean;
    // _filteredMean = _filteredMean + KalmanGain * measurementResidual
    Blas_Mat_Vec_Mult(KalmanGain,measurementResidual,_filteredMean,1.0,1.0);

#ifdef PRINT_INFO
    cout << "KalmanFilter::step: filtered mean after update:" << endl << _filteredMean;
#endif    

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
    // (note that _predictiveCovariance is initalized above to _stateEquationCovariance)
    Blas_Mat_Mat_Trans_Mult(RfilteredCovariance,_R,_predictiveCovariance,1.0,1.0);
}

void KalmanFilter::setFilteredMean(const tVector &filteredMean)
{
    if(filteredMean.size()!=_nElementsToEstimate)    
        throw RuntimeException("KalmanFilter::setFilteredMean: the size of the received vector is wrong.");

    _filteredMean = filteredMean;

    // _predictiveMean = _R*_filteredMean
    Blas_Mat_Vec_Mult(_R,_filteredMean,_predictiveMean);
}

void KalmanFilter::setFilteredCovariance(const tMatrix &filteredCovariance)
{
    if(filteredCovariance.rows()!=_nElementsToEstimate || filteredCovariance.cols()!=_nElementsToEstimate)    
        throw RuntimeException("KalmanFilter::setFilteredCovariance: the dimensions of the received matrix are wrong.");
        
    _filteredCovariance = filteredCovariance;

    tMatrix RfilteredCovariance(_nElementsToEstimate,_nElementsToEstimate);
    
    // RfilteredCovariance = _R*_filteredCovariance
    Blas_Mat_Mat_Mult(_R,_filteredCovariance,RfilteredCovariance);

    _predictiveCovariance = _stateEquationCovariance;

    // _predictiveCovariance = _predictiveCovariance + RfilteredCovariance*_R'
    // (note that _predictiveCovariance is initalized above to _stateEquationCovariance)
    Blas_Mat_Mat_Trans_Mult(RfilteredCovariance,_R,_predictiveCovariance,1.0,1.0); 
}
