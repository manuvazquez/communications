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
_nElementsToEstimate(initialMean.size()),_R(R),_stateEquationCovariance(stateEquationCovariance),_predictiveCovariance(R.rows(),R.rows()),_predictiveMean(R.rows())
{	
	// matrices R and F dimensions are not coherent with those of the vector to be estimated
	assert(R.rows()==_nElementsToEstimate && _nElementsToEstimate==R.cols());
	
	// the number of rows and columns of the covariance must be the number of rows of the mean
	assert(initialMean.size()==initialCovariance.rows() && initialMean.size()==initialCovariance.cols());

    setFilteredMean(initialMean);
    setFilteredCovariance(initialCovariance);
}

void KalmanFilter::step(const MatrixXd &F,const VectorXd &observation,const MatrixXd &observationEquationCovariance)
{
	// the matrix F or observation vector dimensions are wrong
	assert(F.cols()==_nElementsToEstimate && F.rows()==observation.size());
    
    MatrixXd predictiveCovariance_Ft = _predictiveCovariance*F.transpose();
    
    MatrixXd residualCovariance = observationEquationCovariance + F*predictiveCovariance_Ft;

    MatrixXd KalmanGain = predictiveCovariance_Ft*residualCovariance.inverse();
    
    VectorXd measurementResidual = observation - F*_predictiveMean;

    _filteredMean = _predictiveMean + KalmanGain*measurementResidual;

    _filteredCovariance = _predictiveCovariance - KalmanGain*F*_predictiveCovariance;

    // **************** PREDICTIVES **********************
    
    _predictiveMean = _R*_filteredMean;
    _predictiveCovariance = _stateEquationCovariance + _R*_filteredCovariance*_R.transpose();
}

void KalmanFilter::setFilteredMean(const VectorXd &filteredMean)
{
	// the size of the received vector is wrong
	assert(filteredMean.size()==_nElementsToEstimate);

    _filteredMean = filteredMean;

    _predictiveMean = _R*_filteredMean;
}

void KalmanFilter::setFilteredCovariance(const MatrixXd &filteredCovariance)
{
	// the dimensions of the received matrix are wrong
	assert(filteredCovariance.rows()==_nElementsToEstimate && filteredCovariance.cols()==_nElementsToEstimate);
        
    _filteredCovariance = filteredCovariance;
    
    _predictiveCovariance = _stateEquationCovariance + _R*_filteredCovariance*_R.transpose();
}
