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

RLSEstimator::RLSEstimator(tMatrix &initialEstimation,double forgettingFactor): ChannelMatrixEstimator(initialEstimation),_forgettingFactor(forgettingFactor),_invForgettingFactor(1.0/forgettingFactor),_invRtilde(LaGenMatDouble::eye(_Nm)),_pTilde(LaGenMatDouble::zeros(_L,_Nm)),
// auxiliary variables initialization
_invForgettingFactorSymbolsVectorInvRtilde(_Nm),_g(_Nm),_invForgettingFactorInvRtildeSymbolsVector(_Nm),_invForgettingFactorInvRtildeSymbolsVectorg(_Nm,_Nm),_observationsSymbolsVector(_L,_Nm)
// de la otra forma
// ,_transposeEstimatedChannelMatrix(LaGenMatDouble::zeros(_Nm,_L)),_P(LaGenMatDouble::eye(_Nm)),_PsymbolsVector(_Nm),_k(_Nm),_symbolsVectorTransposeEstimatedChannelMatrix(_L),_estimationError(_L),_kEstimationError(_Nm,_L),_kSymbolsVector(_Nm,_Nm),_identityMinuskSymbolsVector(_Nm,_Nm),_auxP(_Nm,_Nm)
{
}

ChannelMatrixEstimator* RLSEstimator::Clone()
{
	// it relies on copy constructor
	return new RLSEstimator(*this);
}

tMatrix RLSEstimator::NextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
	if(observations.size()!=_L || (symbolsMatrix.rows()*symbolsMatrix.cols())!=_Nm)
		throw RuntimeException("Observations vector length or symbols matrix dimensions are wrong.");

	tVector symbolsVector = Util::ToVector(symbolsMatrix,columnwise);

    // _invForgettingFactorSymbolsVectorInvRtilde = symbolsVector'*_invRtilde = _invRtilde'*symbolsVector
    Blas_Mat_Trans_Vec_Mult(_invRtilde,symbolsVector,_invForgettingFactorSymbolsVectorInvRtilde,_invForgettingFactor);

    double auxDenominator = 1.0 + Blas_Dot_Prod(_invForgettingFactorSymbolsVectorInvRtilde,symbolsVector);

    _g = _invForgettingFactorSymbolsVectorInvRtilde;
    _g *= (1.0/auxDenominator);

    // _invForgettingFactorInvRtildeSymbolsVector = _invRtilde*symbolsVector
    Blas_Mat_Vec_Mult(_invRtilde,symbolsVector,_invForgettingFactorInvRtildeSymbolsVector,_invForgettingFactor);

	// _invRtilde = _invForgettingFactorInvRtildeSymbolsVector*_g'
	Util::Mult(_invForgettingFactorInvRtildeSymbolsVector,_g,_invForgettingFactorInvRtildeSymbolsVectorg);

	// _invRtilde = _invForgettingFactor*_invRtilde - _invForgettingFactorInvRtildeSymbolsVectorg
	Util::Add(_invRtilde,_invForgettingFactorInvRtildeSymbolsVectorg,_invRtilde,_invForgettingFactor,-1.0);

	// _observationsSymbolsVector = observations*symbolsVector'
	Util::Mult(observations,symbolsVector,_observationsSymbolsVector);

	// _pTilde = _forgettingFactor*_pTilde + _observationsSymbolsVector
	Util::Add(_pTilde,_observationsSymbolsVector,_pTilde,_forgettingFactor);

	// _pTildeInvRtilde = _pTilde*_invRtilde
	Blas_Mat_Mat_Mult(_pTilde,_invRtilde,_lastEstimatedChannelMatrix);

	return _lastEstimatedChannelMatrix;
}

// tMatrix RLSEstimator::NextMatrix2(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
// {
// 	if(observations.size()!=_L || (symbolsMatrix.rows()*symbolsMatrix.cols())!=_Nm)
// 		throw RuntimeException("Observations vector length or symbols matrix dimensions are wrong.");
// 
// 	tVector symbolsVector = Util::ToVector(symbolsMatrix,columnwise);
// 
// 	//de la otra forma
// 	// _PsymbolsVector = _invForgettingFactor * _P * symbolsVector
// 	Blas_Mat_Vec_Mult(_P,symbolsVector,_PsymbolsVector,_invForgettingFactor);
// 
// 	_k = _PsymbolsVector;
// 
// 	_k *= 1.0/(1.0 + Blas_Dot_Prod(symbolsVector,_PsymbolsVector));
// 
// 	// _symbolsVectorTransposeEstimatedChannelMatrix = symbolsVector' * _transposeEstimatedChannelMatrix' = _transposeEstimatedChannelMatrix * symbolsVector
// 	Blas_Mat_Trans_Vec_Mult(_transposeEstimatedChannelMatrix,symbolsVector,_symbolsVectorTransposeEstimatedChannelMatrix);
// 
// 	// _estimationError = observations - _symbolsVectorTransposeEstimatedChannelMatrix
// 	Util::Add(observations,_symbolsVectorTransposeEstimatedChannelMatrix,_estimationError,1.0,-1.0);
// 
// 	// _kEstimationError = _k * _estimationError'
// 	Util::Mult(_k,_estimationError,_kEstimationError);
// 
// 	// _transposeEstimatedChannelMatrix = _transposeEstimatedChannelMatrix + _kEstimationError
// 	Util::Add(_transposeEstimatedChannelMatrix,_kEstimationError,_transposeEstimatedChannelMatrix);
// 
// 	// _kSymbolsVector = _k * symbolsVector'
// 	Util::Mult(_k,symbolsVector,_kSymbolsVector);
// 
// 	// _identityMinuskSymbolsVector = eye(_Nm) - _kSymbolsVector
// 	Util::Add(LaGenMatDouble::eye(_Nm,_Nm),_kSymbolsVector,_identityMinuskSymbolsVector,1.0,-1.0);
// 
// 	// _auxP = _invForgettingFactor * _identityMinuskSymbolsVector * _P
// 	Blas_Mat_Mat_Mult(_identityMinuskSymbolsVector,_P,_auxP,_invForgettingFactor);
// 
// 	_P = _auxP;
// 
// 	Util::Transpose(_transposeEstimatedChannelMatrix,_lastEstimatedChannelMatrix);
// 
// 	return _lastEstimatedChannelMatrix;
// }
