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
#ifndef DECORRELATORDETECTOR_H
#define DECORRELATORDETECTOR_H

#include <LinearDetector.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>

#include <Eigen/LU>

class DecorrelatorDetector : public LinearDetector
{
protected:
    MatrixXd _filter_eigen;   
public:
    DecorrelatorDetector(int rows, int cols, double alphabetVariance);

    virtual double nthSymbolVariance(int n);
    virtual LinearDetector* clone();
    virtual MatrixXd computedFilter_eigen() { return _filter_eigen;} // eigen
    virtual VectorXd detect(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance);
	virtual void stateStep(VectorXd observations) {}

};

#endif
