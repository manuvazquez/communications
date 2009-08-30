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
#ifndef CDMAKALMANESTIMATOR_H
#define CDMAKALMANESTIMATOR_H

#include <KalmanEstimator.h>

/**
It implements a channel matrix estimator for a Multiuser CDMA autoregressive channel. It assumes a single receiving antenna.

	@author Manu <manu@rustneversleeps>
*/
class CDMAKalmanEstimator : public KalmanEstimator
{
protected:
    MatrixXd _spreadingCodes;
    
    virtual tMatrix buildMeasurementMatrix(const tVector& symbolsVector)
    {
        return Util::eigen2lapack(buildMeasurementMatrix(Util::lapack2eigen(symbolsVector)));
    }
    virtual MatrixXd buildMeasurementMatrix(const VectorXd& symbolsVector); // eigen
public:
    CDMAKalmanEstimator(const MatrixXd& initialEstimation, const MatrixXd& variances, vector< double > ARcoefficients, double ARvariance, const tMatrix &spreadingCodes);

    CDMAKalmanEstimator(const CDMAKalmanEstimator& cdmaKalmanEstimator);
    virtual CDMAKalmanEstimator* clone() const;
    virtual tMatrix sampleFromPredictive() const { return Util::eigen2lapack(sampleFromPredictive_eigen()); }
    virtual MatrixXd sampleFromPredictive_eigen() const; // eigen
};

#endif
