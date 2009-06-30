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

KnownSymbolsKalmanEstimator::KnownSymbolsKalmanEstimator(const tMatrix& initialEstimation, const tMatrix& variances, int N, double ARcoefficient, double ARvariance,const tMatrix &symbols,int startDetectionTime): KalmanEstimator(initialEstimation, variances, N, ARcoefficient, ARvariance),_presentTime(startDetectionTime),_symbols(symbols)
{
}

tMatrix KnownSymbolsKalmanEstimator::nextMatrix(const tVector &observations, const tMatrix &symbolsMatrix, double noiseVariance)
{
	_presentTime++;
	return KalmanEstimator::nextMatrix(observations, _symbols(tRange(0,_N-1),tRange(_presentTime-_m,_presentTime-1)), noiseVariance);
}

KnownSymbolsKalmanEstimator* KnownSymbolsKalmanEstimator::Clone() const
{
	return new KnownSymbolsKalmanEstimator(*this);
}
