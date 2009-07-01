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
#include "CDMAKalmanEstimator.h"

// #define DEBUG

CDMAKalmanEstimator::CDMAKalmanEstimator(tMatrix initialEstimation, const tMatrix &spreadingCodes, const vector<double> &ARprocCoeffs, double ARprocVar): ChannelMatrixEstimator(initialEstimation, spreadingCodes.cols()),_stateVectorLength(initialEstimation.cols()*ARprocCoeffs.size()),_spreadingCodes(spreadingCodes)
{
    if(initialEstimation.rows()!=1)
        throw RuntimeException("CDMAKalmanEstimator::CDMAKalmanEstimator: initial channel estimation is not a row vector.");

    if(initialEstimation.cols()!=_N)
        throw RuntimeException("CDMAKalmanEstimator::CDMAKalmanEstimator: channel is not flat.");

    // _L is computed in the super class as the number of rows in "initialEstimation", which in this case is not correct
    _L = spreadingCodes.rows();

    // stateTransitionMatrix is a blockwise matrix that represents the state transition matrix
    tMatrix stateTransitionMatrix = LaGenMatDouble::zeros(_stateVectorLength,_stateVectorLength);
    
    uint i,j,iBlockCol,iBlockRow;
    
    // state equation matrix is set
    for(iBlockRow=0;iBlockRow<ARprocCoeffs.size()-1;iBlockRow++)
        for(iBlockCol=iBlockRow+1;iBlockCol<ARprocCoeffs.size();iBlockCol++)
            for(j=iBlockCol*initialEstimation.cols(),i=0;i<initialEstimation.cols();j++,i++)
                stateTransitionMatrix(iBlockRow*initialEstimation.cols()+i,j) = 1.0;
    
    for(iBlockCol=0;iBlockCol<ARprocCoeffs.size();iBlockCol++)
        for(j=iBlockCol*initialEstimation.cols(),i=0;i<initialEstimation.cols();j++,i++)
            stateTransitionMatrix((ARprocCoeffs.size()-1)*initialEstimation.cols()+i,j) = ARprocCoeffs[ARprocCoeffs.size()-1-iBlockCol];
    
    // covariance matrix is set
    tMatrix stateEquationCovariance = LaGenMatDouble::zeros(_stateVectorLength,_stateVectorLength);
    
    for(i=_stateVectorLength-initialEstimation.cols();i<_stateVectorLength;i++)
        stateEquationCovariance(i,i) = ARprocVar;

#ifdef DEBUG
    cout << "AR process coefficients ";
    Util::Print(ARprocCoeffs);
    cout << "stateTransitionMatrix = " << endl << stateTransitionMatrix;
    cout << "stateEquationCovariance = " << endl << stateEquationCovariance;
    getchar();
#endif

    tMatrix initialMean = Util::append(initialEstimation,initialEstimation);
    tMatrix initialCovariance = LaGenMatDouble::eye(_stateVectorLength);

    _kalmanFilter = new KalmanFilter(stateTransitionMatrix,stateEquationCovariance,initialMean,initialCovariance);
}


CDMAKalmanEstimator::CDMAKalmanEstimator(int N): ChannelMatrixEstimator(N)
{
}


CDMAKalmanEstimator::~CDMAKalmanEstimator()
{
    delete _kalmanFilter;
}


ChannelMatrixEstimator* CDMAKalmanEstimator::Clone() const
{
}

double CDMAKalmanEstimator::likelihood(const tVector& observations, const tMatrix symbolsMatrix, double noiseVariance)
{
    if(observations.size()!=_L || (symbolsMatrix.rows()*symbolsMatrix.cols())!=_Nm)
        throw RuntimeException("CDMAKalmanEstimator::likelihood: observations or symbols vector length is wrong.");

    // pivots vector needed for factorizations
    tLongIntVector piv(_stateVectorLength);

    tMatrix invPredictiveCovariance = _kalmanFilter->PredictiveCovariance();
    LUFactorizeIP(invPredictiveCovariance,piv);
    
    // detPredictiveCovariance = det(_kalmanFilter->PredictiveCovariance())
    double detPredictiveCovariance = 1.0;
    for(int i=0;i<_stateVectorLength;i++)
        detPredictiveCovariance *= invPredictiveCovariance(i,i);

    // invPredictiveCovariance = inv(_kalmanFilter->PredictiveCovariance())
    LaLUInverseIP(invPredictiveCovariance,piv);

    tMatrix F = BuildFfromSymbolsMatrix(symbolsMatrix);

    tMatrix B = invPredictiveCovariance;
    // B = invPredictiveCovariance + (1.0/noiseVariance) F' * F
    Blas_Mat_Trans_Mat_Mult(F,F,B,1.0/noiseVariance,1.0);


    // B = inv(B) <------------------------------------------------------------------------
    LUFactorizeIP(B,piv);
    LaLUInverseIP(B,piv);

    tVector invPredictiveCovariancePredictiveMean(_stateVectorLength);

    // invPredictiveCovariancePredictiveMean = invPredictiveCovariance * _kalmanFilter->PredictiveMean()
    Blas_Mat_Vec_Mult(invPredictiveCovariance,_kalmanFilter->PredictiveMean(),invPredictiveCovariancePredictiveMean);

    tVector auxAuxArgExp = invPredictiveCovariancePredictiveMean;
    // auxAuxArgExp = invPredictiveCovariancePredictiveMean + (1.0/noiseVariance) F' * observations
    // auxAuxArgExp is "a" in my thesis, appendix D :)
    Blas_Mat_Trans_Vec_Mult(F,observations,auxAuxArgExp,1.0/noiseVariance,1.0);

    tVector auxAuxArgExpInvB(_stateVectorLength);

    // auxAuxArgExpInvB = B' * auxAuxArgExp = auxAuxArgExp * B
    Blas_Mat_Trans_Vec_Mult(B,auxAuxArgExp,auxAuxArgExpInvB);

    // _auxArgExp = auxAuxArgExpInvB . auxAuxArgExp
    double auxArgExp = Blas_Dot_Prod(auxAuxArgExpInvB,auxAuxArgExp);

    tVector observationsNoiseCovariance = observations;
    observationsNoiseCovariance *= noiseVariance;

    // _observationsNoiseCovarianceObservations = observationsNoiseCovariance . observations
    double observationsNoiseCovarianceObservations = Blas_Dot_Prod(observationsNoiseCovariance,observations);

    // predictiveMeanInvPredictiveCovariancePredictiveMean = _kalmanFilter->PredictiveMean() . invPredictiveCovariancePredictiveMean
    double predictiveMeanInvPredictiveCovariancePredictiveMean = Blas_Dot_Prod(_kalmanFilter->PredictiveMean(),invPredictiveCovariancePredictiveMean);

    double argExp = -0.5*(observationsNoiseCovarianceObservations + predictiveMeanInvPredictiveCovariancePredictiveMean - auxArgExp);

    //detInvB = det(B) (recall B = inv(B))
    LUFactorizeIP(B,piv);
    double detInvB = 1.0;
    for(int i=0;i<_stateVectorLength;i++)
        detInvB *= B(i,i);

    return sqrt(fabs(detInvB))/(pow(2*M_PI*noiseVariance,_L/2)*sqrt(fabs(detPredictiveCovariance)))*exp(argExp);    
}

tMatrix CDMAKalmanEstimator::nextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
    if(observations.size()!=_L)
        throw RuntimeException("CDMAKalmanEstimator::nextMatrix: observations vector length is wrong.");

    if(symbolsMatrix.rows()!= _N || symbolsMatrix.cols()!=1)
        throw RuntimeException("CDMAKalmanEstimator::nextMatrix: symbols vector dimensions are wrong.");

    tMatrix observationEquationCovariance = LaGenMatDouble::eye(_L);
    observationEquationCovariance *= noiseVariance;
    _kalmanFilter->Step(BuildFfromSymbolsMatrix(symbolsMatrix),observations,observationEquationCovariance);

    _lastEstimatedChannelMatrix = _kalmanFilter->FilteredMean();

    return  _lastEstimatedChannelMatrix;
}

tMatrix CDMAKalmanEstimator::BuildFfromSymbolsMatrix(const tVector &symbolsVector)
{
    int i,j;
    tMatrix res = LaGenMatDouble::zeros(_L,_N);

    if(symbolsVector.size()!=_N)
        throw RuntimeException("CDMAKalmanEstimator::BuildFfromSymbolsMatrix: the number of elements of the received symbols vector is wrong.");

    // stacks the symbols inside symbolsMatrix to construct F
    for(i=0;i<_L;i++)
        for(j=0;j<_N;j++)
            res(i,j) = symbolsVector(j)*_spreadingCodes(i,j);
    return res;
}
