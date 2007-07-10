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

// #define DEBUG

RLSEstimator::RLSEstimator(int N,double forgettingFactor): ChannelMatrixEstimator(N),_invForgettingFactor(1.0/forgettingFactor)
{
}

RLSEstimator::RLSEstimator(const tMatrix &initialEstimation,int N,double forgettingFactor): ChannelMatrixEstimator(initialEstimation,N),_invForgettingFactor(1.0/forgettingFactor),_invRtilde(LaGenMatDouble::eye(_Nm)),_pTilde(LaGenMatDouble::zeros(_L,_Nm))
{
}

ChannelMatrixEstimator* RLSEstimator::Clone() const
{
	// it relies on copy constructor
	return new RLSEstimator(*this);
}

tMatrix RLSEstimator::NextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
	if(observations.size()!=_L || (symbolsMatrix.rows()*symbolsMatrix.cols())!=_Nm)
		throw RuntimeException("RLSEstimator::NextMatrix: Observations vector length or symbols matrix dimensions are wrong.");

	tVector symbolsVector = Util::ToVector(symbolsMatrix,columnwise);

	tVector invForgettingFactorSymbolsVectorInvRtilde(_Nm);
    // invForgettingFactorSymbolsVectorInvRtilde = symbolsVector'*_invRtilde = _invRtilde'*symbolsVector
    Blas_Mat_Trans_Vec_Mult(_invRtilde,symbolsVector,invForgettingFactorSymbolsVectorInvRtilde,_invForgettingFactor);

    double auxDenominator = 1.0 + Blas_Dot_Prod(invForgettingFactorSymbolsVectorInvRtilde,symbolsVector);

#ifdef DEBUG
	cout << "auxDenominator = " << auxDenominator << endl;
#endif

    tVector g = invForgettingFactorSymbolsVectorInvRtilde;
    g *= (1.0/auxDenominator);

#ifdef DEBUG
	cout << "g" << endl << g;
#endif

	tVector invForgettingFactorInvRtildeSymbolsVector(_Nm);

    // invForgettingFactorInvRtildeSymbolsVector = _invForgettingFactor*_invRtilde*symbolsVector
    Blas_Mat_Vec_Mult(_invRtilde,symbolsVector,invForgettingFactorInvRtildeSymbolsVector,_invForgettingFactor);

    // _invRtilde = _invForgettingFactor*_invRtilde
    _invRtilde *= _invForgettingFactor;

    // _invRtilde = _invRtilde - invForgettingFactorInvRtildeSymbolsVector*g
    Blas_R1_Update(_invRtilde,invForgettingFactorInvRtildeSymbolsVector,g,-1.0);

    // _pTilde = _forgettingFactor*_pTilde
    _pTilde *= (1.0/_invForgettingFactor);

#ifdef DEBUG
	cout << "_pTilde antes" << endl << _pTilde;
	cout << "observations" << endl << observations;
	cout << "symbolsVector" << endl << symbolsVector;
#endif

    // _pTilde = _pTilde + observations*symbolsVector'
    Blas_R1_Update(_pTilde,observations,symbolsVector);

#ifdef DEBUG
	cout << "_pTildeDespues" << endl << _pTilde;
#endif

	// _lastEstimatedChannelMatrix = _pTilde*_invRtilde
	Blas_Mat_Mat_Mult(_pTilde,_invRtilde,_lastEstimatedChannelMatrix);

#ifdef DEBUG
	cout << "Voy a devolver" << endl << _lastEstimatedChannelMatrix;
	cout << "Una tecla..."; getchar();
#endif

	return _lastEstimatedChannelMatrix;
}
