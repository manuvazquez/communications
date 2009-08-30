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

NLMSEstimator::NLMSEstimator(const MatrixXd& initialEstimation, int N, double mu): LMSEstimator(initialEstimation, N, mu)
{
}

NLMSEstimator* NLMSEstimator::clone() const
{
    return new NLMSEstimator(*this);
}

MatrixXd NLMSEstimator::nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance)
{
    VectorXd symbolsVector = Util::toVector(symbolsMatrix,columnwise);
    
    VectorXd error = _lastEstimatedChannelMatrix_eigen*symbolsVector-observations;
    
    _lastEstimatedChannelMatrix_eigen = _lastEstimatedChannelMatrix_eigen - _mu/error.dot(error)*error*symbolsVector.transpose();
    
    return _lastEstimatedChannelMatrix_eigen;
}
