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
#ifndef KNOWNSYMBOLSKALMANESTIMATOR_H
#define KNOWNSYMBOLSKALMANESTIMATOR_H

#include <KalmanEstimator.h>

/**
    @author Manu <manu@rustneversleeps>
*/
class KnownSymbolsKalmanEstimator : public KalmanEstimator
{
protected:
    int _presentTime;
    const MatrixXd &_symbols;
//     const MatrixXd &_symbols_eigen;   
public:
    KnownSymbolsKalmanEstimator(const MatrixXd& initialEstimation, const MatrixXd& variances, int N, vector<double> ARcoefficient, double ARvariance,const MatrixXd &symbols,int startDetectionTime);

    KnownSymbolsKalmanEstimator* clone() const;

//     virtual tMatrix nextMatrix(const tVector &observations, const MatrixXd &symbolsMatrix, double noiseVariance);
    virtual tMatrix nextMatrix(const tVector &observations, const MatrixXd &symbolsMatrix, double noiseVariance)
    {
        return Util::eigen2lapack(nextMatrix(Util::lapack2eigen(observations),symbolsMatrix,noiseVariance));
    }
    virtual MatrixXd nextMatrix(const VectorXd &observations, const MatrixXd &symbolsMatrix, double noiseVariance); // eigen


};

#endif
