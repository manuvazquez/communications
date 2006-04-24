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

RLSEstimator::RLSEstimator(int nRows, int nColumns,double forgettingFactor): ChannelMatrixEstimator(nRows, nColumns),_forgettingFactor(forgettingFactor),_invForgettingFactor(1.0/forgettingFactor),_invRtilde(LaGenMatDouble::eye(_Nm)),_pTilde(LaGenMatDouble::zeros(_L,_Nm)),
// auxiliary variables initialization
_symbolsVectorInvRtilde(_Nm),_g(_Nm),_invRtildeSymbolsVector(_Nm)
{
}

ChannelMatrixEstimator* RLSEstimator::Clone()
{
}

tMatrix RLSEstimator::NextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
	if(observations.size()!=_L || (symbolsMatrix.rows()*symbolsMatrix.cols())!=_Nm)
		throw RuntimeException("Observations vector length or symbols matrix dimensions are wrong.");

	tVector symbolsVector = Util::ToVector(symbolsMatrix,columnwise);

    // _symbolsVectorInvRtilde = symbolsVector'*_invRtilde = _invRtilde'*symbolsVector
    Blas_Mat_Trans_Vec_Mult(_invRtilde,symbolsVector,_symbolsVectorInvRtilde,_invForgettingFactor);

    double auxDenominator = 1.0 + Blas_Dot_Prod(_symbolsVectorInvRtilde,symbolsVector);

    _g = _symbolsVectorInvRtilde;
    _g *= (1.0/auxDenominator);

    // _invRtildeSymbolsVector = _invRtilde*symbolsVector
    Blas_Mat_Vec_Mult(_invRtilde,symbolsVector,_invRtildeSymbolsVector,_invForgettingFactor);

//     Blas_
}

