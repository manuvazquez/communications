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

CDMAKalmanEstimator::CDMAKalmanEstimator(const tMatrix& initialEstimation, const tMatrix& variances, int N, vector< double > ARcoefficients, double ARvariance, const tMatrix &spreadingCodes): KalmanEstimator(initialEstimation, variances, N, ARcoefficients, ARvariance),_spreadingCodes(spreadingCodes)
{
    if(spreadingCodes.cols()!=_nInputs)
        throw RuntimeException("CDMAKalmanEstimator::CDMAKalmanEstimator: the number of spreading codes doesn't match the number of users.");
    
    _nOutputs = _spreadingCodes.rows();
}


CDMAKalmanEstimator::~CDMAKalmanEstimator()
{
}


CDMAKalmanEstimator::CDMAKalmanEstimator(const CDMAKalmanEstimator& cdmaKalmanEstimator):KalmanEstimator(cdmaKalmanEstimator),_spreadingCodes(cdmaKalmanEstimator._spreadingCodes)
{
}

CDMAKalmanEstimator* CDMAKalmanEstimator::Clone() const
{
    return new CDMAKalmanEstimator(*this);
}

tMatrix CDMAKalmanEstimator::BuildFfromSymbolsMatrix(const tVector& symbolsVector)
{
    if(symbolsVector.size()!=_nInputs)
        throw RuntimeException("CDMAKalmanEstimator::BuildFfromSymbolsMatrix: symbols vector length is wrong.");

    tMatrix CS(_nOutputs,_nInputs);
    
    for(uint i=0;i<_nOutputs;i++)
        for(uint j=0;j<_nInputs;j++)
            CS(i,j) = _spreadingCodes(i,j)*symbolsVector(j);
            
    return CS;
}

