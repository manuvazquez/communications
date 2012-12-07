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


#include "AugmentedObservationsKalmanEstimator.h"

AugmentedObservationsKalmanEstimator::AugmentedObservationsKalmanEstimator(const MatrixXd& initialEstimation, const MatrixXd& variances, uint N, std::vector<double> ARcoefficients, double ARvariance):
KalmanEstimator(initialEstimation,variances,N,ARcoefficients,ARvariance),_nARcoeffs(ARcoefficients.size())
// ,_previousSymbolVectors(MatrixXd::Zero(_nInputs,ARcoefficients.size()-1)),_previousObservations(MatrixXd::Zero(_nOutputs,ARcoefficients.size()-1)),_observationEquationCovariances(ARcoefficients.size()-1,MatrixXd::Zero(_nOutputs,_nOutputs))
{
	if(ARcoefficients.size()>1)
	{
		_previousSymbolVectors = MatrixXd::Zero(_nInputs,ARcoefficients.size()-1);
		_previousObservations = MatrixXd::Zero(_nOutputs,ARcoefficients.size()-1);
		_observationEquationCovariances = std::vector<MatrixXd>(ARcoefficients.size()-1,MatrixXd::Zero(_nOutputs,_nOutputs));
	}
}

AugmentedObservationsKalmanEstimator* AugmentedObservationsKalmanEstimator::clone() const
{
    return new AugmentedObservationsKalmanEstimator(*this);
}

MatrixXd AugmentedObservationsKalmanEstimator::nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, const MatrixXd &observationEquationCovariance)
{
	// if there is only one AR coefficient, we can fall back to the super class implementation
	if(_nARcoeffs==1)
		return KalmanEstimator::nextMatrix(observations,symbolsMatrix,observationEquationCovariance);
	
	assert(observations.size()==_nOutputs);
	assert(symbolsMatrix.size()==_nInputsXchannelOrder);

    // a matrix of zeros from which the observation matrix will be built
    MatrixXd extendedObservationMatrix = MatrixXd::Zero(_nOutputs*_nARcoeffs,_nExtStateVectorCoeffs);
	
	// the new observation is stored BELOW the previous ones in a new vector
	VectorXd concatObservations(_previousObservations.rows()+observations.rows(),_previousObservations.cols());
	concatObservations << _previousObservations,observations;
	
	// the symbols vectors are stored TO THE RIGHT of the previous ones in a new matrix
	MatrixXd concatSymbols(_previousSymbolVectors.rows(),_previousSymbolVectors.cols()+symbolsMatrix.cols());
	concatSymbols << _previousSymbolVectors,symbolsMatrix;
	
	// the new observation equation covariance is added to the (c++) vector of previous ones
	_observationEquationCovariances.push_back(observationEquationCovariance);
	
	// the observations matrix is built from the symbol vectors
	for(uint iAR=0;iAR<_nARcoeffs;iAR++)
		extendedObservationMatrix.block(0+iAR*_nOutputs,_nExtStateVectorCoeffs-_nChannelCoeffs*(_nARcoeffs-iAR),_nOutputs,_nChannelCoeffs) = buildObservationMatrix(Util::toVector(concatSymbols.block(0,iAR,_nInputs,symbolsMatrix.cols()),columnwise));
	
	// a KF step is taken using the previously built stuff
	_kalmanFilter->step(extendedObservationMatrix,concatObservations,Util::diag(_observationEquationCovariances));
	
	// the newest observations vectors,...
	_previousObservations = concatObservations.tail((_nARcoeffs-1)*_nOutputs);
	
	// ...symbol vectors...
	_previousSymbolVectors = concatSymbols.block(0,1,_nInputs,_nARcoeffs-1);
	
	// ...and covariance matrices are kept
	_observationEquationCovariances.erase(_observationEquationCovariances.begin());
	
    // notice that only the last coefficients (those representing the channel matrix at current time) are picked up to build the estimated channel matrix
    _lastEstimatedChannelCoefficientsMatrix = Util::toMatrix(_kalmanFilter->filteredMean().tail(_nChannelCoeffs),rowwise,_nChannelMatrixRows);
	
    return _lastEstimatedChannelCoefficientsMatrix;
}

double AugmentedObservationsKalmanEstimator::likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance)
{
    throw RuntimeException("AugmentedObservationsKalmanEstimator::likelihood: not implemented!!");
}

