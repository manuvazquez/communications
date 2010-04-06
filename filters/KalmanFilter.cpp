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
// #define DEBUG
 
KalmanFilter::KalmanFilter(const MatrixXd &R,const MatrixXd &stateEquationCovariance,const VectorXd &initialMean,const MatrixXd &initialCovariance):
_nElementsToEstimate(initialMean.size()),_R_eigen(R),_stateEquationCovariance_eigen(stateEquationCovariance),_predictiveCovariance_eigen(R.rows(),R.rows()),_predictiveMean_eigen(R.rows())
{
    if(R.rows()!=_nElementsToEstimate || _nElementsToEstimate!=R.cols())
        throw RuntimeException("KalmanFilter::KalmanFilter: matrices R and F dimensions are not coherent with those of the vector to be estimated.");

    if(initialMean.size()!=initialCovariance.rows() || initialMean.size()!=initialCovariance.cols())
        throw RuntimeException("KalmanFilter::KalmanFilter: the number of rows and columns of the covariance must be the number of rows of the mean.");

    setFilteredMean(initialMean);
    setFilteredCovariance(initialCovariance);
}

void KalmanFilter::step(const MatrixXd &F_eigen,const VectorXd &observation_eigen,const MatrixXd &observationEquationCovariance_eigen)
{
    if(F_eigen.cols()!=_nElementsToEstimate || F_eigen.rows()!=observation_eigen.size())
        throw RuntimeException("The matrix F or observation vector dimensions are wrong.");
    
    MatrixXd predictiveCovariance_Ft_eigen = _predictiveCovariance_eigen*F_eigen.transpose();
    
    MatrixXd residualCovariance_eigen = observationEquationCovariance_eigen + F_eigen*predictiveCovariance_Ft_eigen;

    MatrixXd KalmanGain_eigen = predictiveCovariance_Ft_eigen*residualCovariance_eigen.inverse();
    
    VectorXd measurementResidual_eigen = observation_eigen - F_eigen*_predictiveMean_eigen;

    _filteredMean_eigen = _predictiveMean_eigen + KalmanGain_eigen*measurementResidual_eigen;

    _filteredCovariance_eigen = _predictiveCovariance_eigen - KalmanGain_eigen*F_eigen*_predictiveCovariance_eigen;

    // **************** PREDICTIVES **********************
    
    _predictiveMean_eigen = _R_eigen*_filteredMean_eigen;
    _predictiveCovariance_eigen = _stateEquationCovariance_eigen + _R_eigen*_filteredCovariance_eigen*_R_eigen.transpose();
}

void KalmanFilter::setFilteredMean(const VectorXd &filteredMean)
{
    if(filteredMean.size()!=_nElementsToEstimate)    
        throw RuntimeException("KalmanFilter::setFilteredMean: the size of the received vector is wrong.");

    _filteredMean_eigen = filteredMean;

    _predictiveMean_eigen = _R_eigen*_filteredMean_eigen;
}

void KalmanFilter::setFilteredCovariance(const MatrixXd &filteredCovariance)
{
    if(filteredCovariance.rows()!=_nElementsToEstimate || filteredCovariance.cols()!=_nElementsToEstimate)    
        throw RuntimeException("KalmanFilter::setFilteredCovariance: the dimensions of the received matrix are wrong.");
        
    _filteredCovariance_eigen = filteredCovariance;
    
    _predictiveCovariance_eigen = _stateEquationCovariance_eigen + _R_eigen*_filteredCovariance_eigen*_R_eigen.transpose();
}
