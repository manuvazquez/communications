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

RLSEstimator::RLSEstimator(const tMatrix &initialEstimation,int N,double forgettingFactor): ChannelMatrixEstimator(initialEstimation,N),_invForgettingFactor(1.0/forgettingFactor),_invRtilde_eigen(MatrixXd::Identity(_nInputsXchannelOrder,_nInputsXchannelOrder))
{
}

ChannelMatrixEstimator* RLSEstimator::clone() const
{
	// it relies on copy constructor
	return new RLSEstimator(*this);
}

// eigen
double RLSEstimator::likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance)
{
    return StatUtil::normalPdf(observations,_lastEstimatedChannelMatrix_eigen*Util::toVector(symbolsMatrix,columnwise),noiseVariance);
}

// eigen
MatrixXd RLSEstimator::nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance)
{
    if(observations.size()!=_nOutputs || symbolsMatrix.size()!=_nInputsXchannelOrder)
        throw RuntimeException("RLSEstimator::NextMatrix: Observations vector length or symbols matrix dimensions are wrong.");

    VectorXd symbolsVector = Util::toVector(symbolsMatrix,columnwise);

    VectorXd invForgettingFactorSymbolsVectorInvRtilde = _invForgettingFactor*_invRtilde_eigen.transpose()*symbolsVector;

    VectorXd g = invForgettingFactorSymbolsVectorInvRtilde / (1.0 + symbolsVector.dot(invForgettingFactorSymbolsVectorInvRtilde));

    _lastEstimatedChannelMatrix_eigen = _lastEstimatedChannelMatrix_eigen + (observations-_lastEstimatedChannelMatrix_eigen*symbolsVector)*g.transpose();

    _invRtilde_eigen = _invForgettingFactor*_invRtilde_eigen - (_invForgettingFactor*_invRtilde_eigen*symbolsVector)*g.transpose();

    _lastEstimatedChannelMatrix = Util::eigen2lapack(_lastEstimatedChannelMatrix_eigen);
    return _lastEstimatedChannelMatrix_eigen;
}
