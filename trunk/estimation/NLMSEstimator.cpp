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
#include "NLMSEstimator.h"

NLMSEstimator::NLMSEstimator(const tMatrix& initialEstimation, int N, double mu): LMSEstimator(initialEstimation, N, mu)
{
}

LMSEstimator* NLMSEstimator::Clone() const
{
    return new NLMSEstimator(*this);
}

tMatrix NLMSEstimator::NextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
    tVector symbolsVector = Util::ToVector(symbolsMatrix,columnwise);

    // _error = observations
    tVector error = observations;

    // error = _lastEstimatedChannelMatrix*symbolsVector - error
    // (note that error is initialized to observations, so that:
    // error = _lastEstimatedChannelMatrix*symbolsVector - observations)
    Blas_Mat_Vec_Mult(_lastEstimatedChannelMatrix,symbolsVector,error,1.0,-1.0);

    // _lastEstimatedChannelMatrix = - _mu * _error * symbolsVector' + _lastEstimatedChannelMatrix
    Blas_R1_Update(_lastEstimatedChannelMatrix,error,symbolsVector,-_mu/Blas_Dot_Prod(error,error));

    return _lastEstimatedChannelMatrix;
}
