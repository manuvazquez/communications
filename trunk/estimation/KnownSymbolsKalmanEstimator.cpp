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
#include "KnownSymbolsKalmanEstimator.h"

// #define DEBUG

KnownSymbolsKalmanEstimator::KnownSymbolsKalmanEstimator(const tMatrix& initialEstimation, const tMatrix& variances, int N, vector<double> ARcoefficients, double ARvariance,const tMatrix &symbols,int startDetectionTime): KalmanEstimator(initialEstimation, variances, N, ARcoefficients, ARvariance),_presentTime(startDetectionTime),_symbols(symbols)
// ,_symbols_eigen(Util::lapack2eigen(symbols))
{
}

// eigen
MatrixXd KnownSymbolsKalmanEstimator::nextMatrix(const VectorXd &observations, const MatrixXd &symbolsMatrix, double noiseVariance)
{
    _presentTime++;
    return KalmanEstimator::nextMatrix(observations, Util::lapack2eigen(_symbols).block(0,_presentTime-_channelOrder,_nInputs,_channelOrder), noiseVariance);
//     return KalmanEstimator::nextMatrix(observations, _symbols(tRange(0,_nInputs-1),tRange(_presentTime-_channelOrder,_presentTime-1)), noiseVariance);
}

KnownSymbolsKalmanEstimator* KnownSymbolsKalmanEstimator::clone() const
{
	return new KnownSymbolsKalmanEstimator(*this);
}
