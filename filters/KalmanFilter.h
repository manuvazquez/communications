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
#ifndef KALMANFILTER_H
#define KALMANFILTER_H

/**
    @author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <exceptions.h>
#include <Util.h>

#include <Eigen/Core>
#include <Eigen/LU>

class KalmanFilter{
private:
    uint _nElementsToEstimate;
    MatrixXd _R,_stateEquationCovariance,_predictiveCovariance,_filteredCovariance;
    VectorXd _predictiveMean,_filteredMean;

public:
    KalmanFilter(const MatrixXd &R,const MatrixXd &stateEquationCovariance,const VectorXd &initialMean,const MatrixXd &initialCovariance);

    void step(const MatrixXd &F,const VectorXd &observation,const MatrixXd &observationEquationCovariance);
    
    VectorXd predictiveMean() const { return _predictiveMean;}
    VectorXd filteredMean() const { return _filteredMean;}
    MatrixXd predictiveCovariance() const { return _predictiveCovariance;}
    MatrixXd filteredCovariance() const { return _filteredCovariance;}
    
    void setFilteredMean(const VectorXd &filteredMean);
    void setFilteredCovariance(const MatrixXd &filteredCovariance);    

};

#endif
