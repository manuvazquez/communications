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
#include "RLSEstimator.h"

RLSEstimator::RLSEstimator(const tMatrix &initialEstimation,int N,double forgettingFactor): ChannelMatrixEstimator(initialEstimation,N),_invForgettingFactor(1.0/forgettingFactor),_invRtilde(LaGenMatDouble::eye(_nInputsXchannelOrder))
{
}

ChannelMatrixEstimator* RLSEstimator::Clone() const
{
	// it relies on copy constructor
	return new RLSEstimator(*this);
}

double RLSEstimator::likelihood(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance)
{
	tVector computedObservations(observations.size());
	Blas_Mat_Vec_Mult(_lastEstimatedChannelMatrix,Util::toVector(symbolsMatrix,columnwise),computedObservations);

	return StatUtil::NormalPdf(observations,computedObservations,noiseVariance);
}

tMatrix RLSEstimator::nextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
    if(observations.size()!=_nOutputs || (symbolsMatrix.rows()*symbolsMatrix.cols())!=_nInputsXchannelOrder)
        throw RuntimeException("RLSEstimator::NextMatrix: Observations vector length or symbols matrix dimensions are wrong.");

    tVector symbolsVector = Util::toVector(symbolsMatrix,columnwise);

    tVector invForgettingFactorSymbolsVectorInvRtilde(_nInputsXchannelOrder);
    // invForgettingFactorSymbolsVectorInvRtilde = symbolsVector'*_invRtilde = _invRtilde'*symbolsVector
    Blas_Mat_Trans_Vec_Mult(_invRtilde,symbolsVector,invForgettingFactorSymbolsVectorInvRtilde,_invForgettingFactor);

    double auxDenominator = 1.0 + Blas_Dot_Prod(invForgettingFactorSymbolsVectorInvRtilde,symbolsVector);

    tVector g = invForgettingFactorSymbolsVectorInvRtilde;
    g *= (1.0/auxDenominator);


    tVector observationsMinusPredictedObservations = observations;
    // observationsMinusPredictedObservations = observationsMinusPredictedObservations - _lastEstimatedChannelMatrix * symbolsVector
    Blas_Mat_Vec_Mult(_lastEstimatedChannelMatrix,symbolsVector,observationsMinusPredictedObservations,-1.0,1.0);

    // _lastEstimatedChannelMatrix = _lastEstimatedChannelMatrix + observationsMinusPredictedObservations*g
    Blas_R1_Update(_lastEstimatedChannelMatrix,observationsMinusPredictedObservations,g);

    tVector invForgettingFactorInvRtildeSymbolsVector(_nInputsXchannelOrder);

    // invForgettingFactorInvRtildeSymbolsVector = _invForgettingFactor*_invRtilde*symbolsVector
    Blas_Mat_Vec_Mult(_invRtilde,symbolsVector,invForgettingFactorInvRtildeSymbolsVector,_invForgettingFactor);

    // _invRtilde = _invForgettingFactor*_invRtilde
    _invRtilde *= _invForgettingFactor;

    // _invRtilde = _invRtilde - invForgettingFactorInvRtildeSymbolsVector*g
    Blas_R1_Update(_invRtilde,invForgettingFactorInvRtildeSymbolsVector,g,-1.0);

    return _lastEstimatedChannelMatrix;
}

