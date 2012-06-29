/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2012  Manu <manuavazquez@gmail.com>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#ifndef AUGMENTEDOBSERVATIONSKALMANESTIMATOR_H
#define AUGMENTEDOBSERVATIONSKALMANESTIMATOR_H

#include <KalmanEstimator.h>


class AugmentedObservationsKalmanEstimator : public KalmanEstimator
{
protected:
	const uint _nARcoeffs;
	MatrixXd _previousSymbolVectors;
	VectorXd _previousObservations;
	std::vector<MatrixXd> _observationEquationCovariances;
	
//     virtual MatrixXd buildObservationMatrix(const VectorXd& symbolsVector);

public:
    AugmentedObservationsKalmanEstimator(const MatrixXd& initialEstimation, const MatrixXd& variances, uint N, std::vector<double> ARcoefficients, double ARvariance);
    virtual AugmentedObservationsKalmanEstimator* clone() const;
    virtual MatrixXd nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, const MatrixXd &observationEquationCovariance);
    virtual double likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance);
};

#endif // AUGMENTEDOBSERVATIONSKALMANESTIMATOR_H
