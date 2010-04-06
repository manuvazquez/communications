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
#include "DecorrelatorDetector.h"

DecorrelatorDetector::DecorrelatorDetector(int rows, int cols, double alphabetVariance): LinearDetector(rows, cols, alphabetVariance)
{
}

double DecorrelatorDetector::nthSymbolVariance(int n,double noiseVariance)
{
//     return _filter_eigen.row(n).dot(_filter_eigen.row(n));
    return noiseVariance*_filter_eigen.row(n).dot(_filter_eigen.row(n));    
}

LinearDetector* DecorrelatorDetector::clone()
{
    return new DecorrelatorDetector(*this);
}

VectorXd DecorrelatorDetector::detect(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance)
{
    _filter_eigen = (channelMatrix.transpose()*channelMatrix).inverse()*channelMatrix.transpose();

    return _filter_eigen*observations;
}
