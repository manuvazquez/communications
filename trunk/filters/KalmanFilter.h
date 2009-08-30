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
#include <lapackpp/gmd.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>

#include <Eigen/Core>
#include <Eigen/LU>
// USING_PART_OF_NAMESPACE_EIGEN

class KalmanFilter{
private:
    int _nElementsToEstimate;
    MatrixXd _R_eigen,_stateEquationCovariance_eigen,_predictiveCovariance_eigen,_filteredCovariance_eigen;
    VectorXd _predictiveMean_eigen,_filteredMean_eigen;

public:
    KalmanFilter(const MatrixXd &R,const MatrixXd &stateEquationCovariance,const VectorXd &initialMean,const MatrixXd &initialCovariance); // eigen

    void step(const MatrixXd &F_eigen,const VectorXd &observation_eigen,const MatrixXd &observationEquationCovariance_eigen); // eigen

    tVector predictiveMean() const { return Util::eigen2lapack(_predictiveMean_eigen);}
    tVector filteredMean() const { return Util::eigen2lapack(_filteredMean_eigen);}
    tMatrix predictiveCovariance() const { return Util::eigen2lapack(_predictiveCovariance_eigen);}
    tMatrix filteredCovariance() const { return Util::eigen2lapack(_filteredCovariance_eigen);}
    
    VectorXd predictiveMean_eigen() const { return _predictiveMean_eigen;} // eigen
    VectorXd filteredMean_eigen() const { return _filteredMean_eigen;} // eigen
    MatrixXd predictiveCovariance_eigen() const { return _predictiveCovariance_eigen;} // eigen
    MatrixXd filteredCovariance_eigen() const { return _filteredCovariance_eigen;} //eigen
    
    void setFilteredMean(const VectorXd &filteredMean);
    void setFilteredCovariance(const MatrixXd &filteredCovariance);    

};

#endif
