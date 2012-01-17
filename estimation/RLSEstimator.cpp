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

RLSEstimator::RLSEstimator(const MatrixXd &initialEstimation,uint N,double forgettingFactor): ChannelMatrixEstimator(initialEstimation,N),_invForgettingFactor(1.0/forgettingFactor),_invRtilde(MatrixXd::Identity(_nInputsXchannelOrder,_nInputsXchannelOrder))
{
}

ChannelMatrixEstimator* RLSEstimator::clone() const
{
	// it relies on copy constructor
	return new RLSEstimator(*this);
}

double RLSEstimator::likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance)
{
//     return StatUtil::normalPdf(observations,_lastEstimatedChannelCoefficientsMatrix*Util::toVector(symbolsMatrix,columnwise),noiseVariance);
	return StatUtil::normalPdf(observations,lastEstimatedChannelMatrix()*Util::toVector(symbolsMatrix,columnwise),noiseVariance);
}

MatrixXd RLSEstimator::nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance)
{
	assert(observations.size()==_nOutputs && symbolsMatrix.size()==_nInputsXchannelOrder);

    VectorXd symbolsVector = Util::toVector(symbolsMatrix,columnwise);

    VectorXd invForgettingFactorSymbolsVectorInvRtilde = _invForgettingFactor*_invRtilde.transpose()*symbolsVector;

    VectorXd g = invForgettingFactorSymbolsVectorInvRtilde / (1.0 + symbolsVector.dot(invForgettingFactorSymbolsVectorInvRtilde));

    _lastEstimatedChannelCoefficientsMatrix = _lastEstimatedChannelCoefficientsMatrix + (observations-_lastEstimatedChannelCoefficientsMatrix*symbolsVector)*g.transpose();

    _invRtilde = _invForgettingFactor*_invRtilde - (_invForgettingFactor*_invRtilde*symbolsVector)*g.transpose();

    return _lastEstimatedChannelCoefficientsMatrix;
}
