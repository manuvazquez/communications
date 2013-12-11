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
    uint _presentTime;
    const MatrixXd &_symbols;
public:
    KnownSymbolsKalmanEstimator(const MatrixXd& initialEstimation, const MatrixXd& variances, uint N, vector<double> ARcoefficient, double ARvariance,const MatrixXd &symbols,uint startDetectionTime);

    KnownSymbolsKalmanEstimator* clone() const;

    virtual MatrixXd nextMatrix(const VectorXd &observations, const MatrixXd &symbolsMatrix, double noiseVariance);


};

#endif
