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
#include "LMSEstimator.h"

LMSEstimator::LMSEstimator(const tMatrix &initialEstimation,int N,double mu): ChannelMatrixEstimator(initialEstimation,N),_mu(mu),_predictedObservations(_L),_error(_L),_deltaMatrix(_L,_Nm)
{
}

LMSEstimator* LMSEstimator::Clone()
{
	return new LMSEstimator(*this);
}

tMatrix LMSEstimator::NextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
	tVector symbolsVector = Util::ToVector(symbolsMatrix,columnwise);

	// _predictedObservations = _lastEstimatedChannelMatrix * symbolsVector
	Blas_Mat_Vec_Mult(_lastEstimatedChannelMatrix,symbolsVector,_predictedObservations);

	// _error = _predictedObservations - observations
	Util::Add(_predictedObservations,observations,_error,1.0,-1.0);

	// _deltaMatrix = _mu * _error * symbolsVector'
	Util::Mult(_error,symbolsVector,_deltaMatrix,_mu);

	// _lastEstimatedChannelMatrix = _lastEstimatedChannelMatrix - _deltaMatrix
	Util::Add(_lastEstimatedChannelMatrix,_deltaMatrix,_lastEstimatedChannelMatrix,1.0,-1.0);

	return _lastEstimatedChannelMatrix;
}

