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

// MatrixXd AugmentedObservationsKalmanEstimator::buildObservationMatrix(const VectorXd& symbolsVector)
// {
//     return KalmanEstimator::buildObservationMatrix(symbolsVector);
// }

AugmentedObservationsKalmanEstimator::AugmentedObservationsKalmanEstimator(const MatrixXd& initialEstimation, const MatrixXd& variances, uint N, std::vector<double> ARcoefficients, double ARvariance):
KalmanEstimator(initialEstimation,variances,N,ARcoefficients,ARvariance),_nARcoeffs(ARcoefficients.size())
// ,_previousSymbolVectors(MatrixXd::Zero(_nInputs,ARcoefficients.size()-1)),_previousObservations(MatrixXd::Zero(_nOutputs,ARcoefficients.size()-1))
// ,_observationEquationCovariances(ARcoefficients.size()-1,MatrixXd::Zero(_nOutputs,_nOutputs))
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
	if(_nARcoeffs==1)
		return KalmanEstimator::nextMatrix(observations,symbolsMatrix,observationEquationCovariance);
	
	assert(observations.size()==_nOutputs);
	assert(symbolsMatrix.size()==_nInputsXchannelOrder);

    // extStateMeasurementMatrix is a matrix of zeros whose right side is the common observation matrix (it is meant to take into account when there is more than one AR coefficient)
    MatrixXd extendedObservationMatrix = MatrixXd::Zero(_nOutputs*_nARcoeffs,_nExtStateVectorCoeffs);
	
	VectorXd concatObservations(_previousObservations.rows()+observations.rows(),_previousObservations.cols());
	concatObservations << _previousObservations,observations;
	
	MatrixXd concatSymbols(_previousSymbolVectors.rows(),_previousSymbolVectors.cols()+symbolsMatrix.cols());
	concatSymbols << _previousSymbolVectors,symbolsMatrix;
	
	_observationEquationCovariances.push_back(observationEquationCovariance);
	
// 	VectorXd concatObservations(_nOutputs*_nARcoeffs);
// 	for(uint iAR=0;iAR<_nARcoeffs;iAR++)
// 		concatObservations.segment(iAR*_nOutputs,_nOutputs) = 
	
	for(uint iAR=0;iAR<_nARcoeffs;iAR++)
		extendedObservationMatrix.block(0+iAR*_nOutputs,_nExtStateVectorCoeffs-_nChannelCoeffs*(_nARcoeffs-iAR),_nOutputs,_nChannelCoeffs) = buildObservationMatrix(Util::toVector(concatSymbols.block(0,iAR,_nInputs,symbolsMatrix.cols()),columnwise));
	
//     extendedObservationMatrix.block(0,_nExtStateVectorCoeffs-_nChannelCoeffs,_nOutputs,_nChannelCoeffs) = buildObservationMatrix(Util::toVector(symbolsMatrix,columnwise));
    
//     _kalmanFilter->step(extendedObservationMatrix,observations,observationEquationCovariance);
		
// 	cout << "extendedObservationMatrix:" << endl << extendedObservationMatrix << endl;
// 	cout << "concatObservations" << endl << concatObservations << endl;
// 	getchar();
	
	_kalmanFilter->step(extendedObservationMatrix,concatObservations,Util::diag(_observationEquationCovariances));
	
	
	_previousObservations = concatObservations.tail((_nARcoeffs-1)*_nOutputs);
	_previousSymbolVectors = concatSymbols.block(0,1,_nInputs,_nARcoeffs-1);
	_observationEquationCovariances.erase(_observationEquationCovariances.begin());
	
    // notice that only the last coefficients (those representing the channel matrix at current time) are picked up to build the estimated channel matrix
    _lastEstimatedChannelCoefficientsMatrix = Util::toMatrix(_kalmanFilter->filteredMean().tail(_nChannelCoeffs),rowwise,_nChannelMatrixRows);
	
    return _lastEstimatedChannelCoefficientsMatrix;
}

double AugmentedObservationsKalmanEstimator::likelihood(const VectorXd& observations, const MatrixXd symbolsMatrix, double noiseVariance)
{
    throw RuntimeException("AugmentedObservationsKalmanEstimator::likelihood: not implemented!!");
}

