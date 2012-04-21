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
#include "KalmanEstimator.h"

// #define PRINT_INFO

KalmanEstimator::KalmanEstimator(const MatrixXd &initialEstimation,const MatrixXd &variances,uint N,vector<double> ARcoefficients,double ARvariance): ChannelMatrixEstimator(initialEstimation,N),_nExtStateVectorCoeffs(_nChannelCoeffs*ARcoefficients.size())
{
    // stateTransitionMatrix is a blockwise matrix that represents the state transition matrix
    MatrixXd stateTransitionMatrix = MatrixXd::Zero(_nExtStateVectorCoeffs,_nExtStateVectorCoeffs);    
    
    uint i,j,iBlockCol,iBlockRow;
    
    // state equation matrix is set
    for(iBlockRow=0;iBlockRow<ARcoefficients.size()-1;iBlockRow++)
        for(iBlockCol=iBlockRow+1;iBlockCol<ARcoefficients.size();iBlockCol++)
            for(j=iBlockCol*_nChannelCoeffs,i=0;i<_nChannelCoeffs;j++,i++)
                stateTransitionMatrix(iBlockRow*_nChannelCoeffs+i,j) = 1.0;
    
    for(iBlockCol=0;iBlockCol<ARcoefficients.size();iBlockCol++)
        for(j=iBlockCol*_nChannelCoeffs,i=0;i<_nChannelCoeffs;j++,i++)
            stateTransitionMatrix((ARcoefficients.size()-1)*_nChannelCoeffs+i,j) = ARcoefficients[ARcoefficients.size()-1-iBlockCol];
    
    // covariance matrix is set
    MatrixXd stateEquationCovariance = MatrixXd::Zero(_nExtStateVectorCoeffs,_nExtStateVectorCoeffs);
    
    for(i=_nExtStateVectorCoeffs-_nChannelCoeffs;i<_nExtStateVectorCoeffs;i++)
        stateEquationCovariance(i,i) = ARvariance;

    VectorXd initialMeanVector(_nExtStateVectorCoeffs);
    MatrixXd initialCovariance = MatrixXd::Zero(_nExtStateVectorCoeffs,_nExtStateVectorCoeffs);
    for(iBlockRow=0;iBlockRow<ARcoefficients.size();iBlockRow++)
        for(i=0;i<_nChannelCoeffs;i++)
        {
            initialMeanVector(iBlockRow*_nChannelCoeffs+i) = initialEstimation(i/_nInputsXchannelOrder,i%_nInputsXchannelOrder);
            initialCovariance(iBlockRow*_nChannelCoeffs+i,iBlockRow*_nChannelCoeffs+i) = variances(i/_nInputsXchannelOrder,i%_nInputsXchannelOrder);
        }

#ifdef PRINT_INFO
    cout << "KalmanEstimator::KalmanEstimator: constructed a Kalman Filter with parameters:" << endl;
    cout << "state transition matrix:" << endl << stateTransitionMatrix;
    cout << "state equation covariance:" << endl << stateEquationCovariance;
    cout << "initial mean" << endl << initialMeanVector;
    cout << "initial covariance" << endl << initialCovariance;
#endif

    _kalmanFilter = new KalmanFilter(stateTransitionMatrix,stateEquationCovariance,initialMeanVector,initialCovariance);    
}

KalmanEstimator::KalmanEstimator(const KalmanEstimator &kalmanEstimator):ChannelMatrixEstimator(kalmanEstimator),_kalmanFilter(new KalmanFilter(*(kalmanEstimator._kalmanFilter))),_nExtStateVectorCoeffs(kalmanEstimator._nExtStateVectorCoeffs)
{
}

KalmanEstimator::~KalmanEstimator()
{
    delete _kalmanFilter;
}

MatrixXd KalmanEstimator::nextMatrix(const VectorXd &observations,const MatrixXd &symbolsMatrix,double noiseVariance)
{
    if(observations.size()!=_nOutputs || symbolsMatrix.size()!=_nInputsXchannelOrder)
	{
	  cout << "observations.size() = " << observations.size() << " _nOutputs = " << _nOutputs << endl;
	  cout << " symbolsMatrix.size() = " << symbolsMatrix.size() << " _nInputsXchannelOrder = " << _nInputsXchannelOrder << endl;
	  throw RuntimeException("KalmanEstimator::NextMatrix: observations vector length or symbols matrix length are wrong.");
	}

    MatrixXd observationEquationCovariance = noiseVariance*MatrixXd::Identity(_nOutputs,_nOutputs);
    
    // extStateMeasurementMatrix is a matrix of zeros whose right side is the common observation matrix. It is meant to take into account when there is more
    // than one AR coefficient
    MatrixXd extStateMeasurementMatrix = MatrixXd::Zero(_nOutputs,_nExtStateVectorCoeffs);    
    
    extStateMeasurementMatrix.block(0,_nExtStateVectorCoeffs-_nChannelCoeffs,_nOutputs,_nChannelCoeffs) = buildMeasurementMatrix(Util::toVector(symbolsMatrix,columnwise));    
    
    _kalmanFilter->step(extStateMeasurementMatrix,observations,observationEquationCovariance);
    
    // notice that only the last coefficients (those representing the channel matrix at current time) are picked up to build the estimated channel matrix
    _lastEstimatedChannelCoefficientsMatrix = Util::toMatrix(_kalmanFilter->filteredMean().tail(_nChannelCoeffs),rowwise,_nChannelMatrixRows);

    return _lastEstimatedChannelCoefficientsMatrix;
}

MatrixXd KalmanEstimator::buildMeasurementMatrix(const VectorXd &symbolsVector)
{
    uint i,j;
    MatrixXd res = MatrixXd::Zero(_nOutputs,_nChannelCoeffs);

    if(symbolsVector.size()!=_nInputsXchannelOrder)
        throw RuntimeException("KalmanEstimator::buildMeasurementMatrix: the number of elements of the received symbols vector is wrong.");

    // stacks the symbols inside symbolsMatrix to construct F
    for(i=0;i<_nOutputs;i++)
        for(j=0;j<_nInputsXchannelOrder;j++)
            res(i,i*_nInputsXchannelOrder+j) = symbolsVector(j);
    return res;
}

double KalmanEstimator::likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance)
{
    if(observations.size()!=_nOutputs || symbolsMatrix.size()!=_nInputsXchannelOrder)
        throw RuntimeException("KalmanEstimator::likelihood: observations vector length or symbols matrix length are wrong.");
        
    Eigen::LDLT<MatrixXd> ldltOfPredictiveCovariance(_kalmanFilter->predictiveCovariance());

    MatrixXd invPredictiveCovariance = MatrixXd::Identity(_nExtStateVectorCoeffs,_nExtStateVectorCoeffs);
    ldltOfPredictiveCovariance.solveInPlace(invPredictiveCovariance);
    
    
    double invPredictiveCovarianceDeterminant = 1.0;    
    for(uint i=0;i<ldltOfPredictiveCovariance.vectorD().rows();i++)
        invPredictiveCovarianceDeterminant *= ldltOfPredictiveCovariance.vectorD().coeff(i);        
        
        
    MatrixXd extStateMeasurementMatrix = MatrixXd::Zero(_nOutputs,_nExtStateVectorCoeffs);
	extStateMeasurementMatrix.block(0,_nExtStateVectorCoeffs-_nChannelCoeffs,_nOutputs,_nChannelCoeffs) = buildMeasurementMatrix(Util::toVector(symbolsMatrix,columnwise));
    
    
	MatrixXd invNoiseVariance_extStateMeasurementMatrixT = 1.0/noiseVariance*extStateMeasurementMatrix.transpose();

	// the LU decomposition is computed...
	PartialPivLU<MatrixXd> luforB(invPredictiveCovariance + invNoiseVariance_extStateMeasurementMatrixT*extStateMeasurementMatrix);

	//...to invert the matrix
	// TODO: is this better than just calling .inv()?
	MatrixXd invB = luforB.solve(MatrixXd::Identity(_nExtStateVectorCoeffs,_nExtStateVectorCoeffs));
    
	VectorXd invPredictiveCovariancePredictiveMean = invPredictiveCovariance*_kalmanFilter->predictiveMean();
    
	VectorXd auxAuxArgExp = invPredictiveCovariancePredictiveMean + invNoiseVariance_extStateMeasurementMatrixT*observations;

	VectorXd auxAuxArgExpInvB = invB.transpose()*auxAuxArgExp;
    
    double argExp = -0.5*(1.0/noiseVariance*observations.dot(observations) + _kalmanFilter->predictiveMean().dot(invPredictiveCovariancePredictiveMean) - auxAuxArgExpInvB.dot(auxAuxArgExp));

    return sqrt(fabs(1.0/luforB.determinant()))/(pow(2*M_PI*noiseVariance,_nOutputs/2)*sqrt(fabs(invPredictiveCovarianceDeterminant)))*exp(argExp);
}

KalmanEstimator *KalmanEstimator::clone() const
{
    // it relies on copy constructor
    return new KalmanEstimator(*this);
}

MatrixXd KalmanEstimator::samplePredicted() const
{
    return Util::toMatrix(
        StatUtil::randnMatrix(_kalmanFilter->predictiveMean(),_kalmanFilter->predictiveCovariance()).tail(_nChannelCoeffs)
        ,rowwise,_nChannelMatrixRows);
}

void KalmanEstimator::setFirstEstimatedChannelMatrix(const MatrixXd &matrix)
{
    ChannelMatrixEstimator::setFirstEstimatedChannelMatrix(matrix);

    VectorXd extState(_nExtStateVectorCoeffs);
    
    uint i=0;
    
    while(i<_nExtStateVectorCoeffs)
    {
        extState(i) = matrix((i%_nChannelCoeffs)/_nInputsXchannelOrder,(i%_nChannelCoeffs)%_nInputsXchannelOrder);
        i++;
    }
     
    _kalmanFilter->setFilteredMean(extState);
}
