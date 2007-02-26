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
#include "RMMSEDetector.h"

RMMSEDetector::RMMSEDetector(int rows, int cols,double alphabetVariance,double forgettingFactor,int nSymbolsToBeDetected): LinearDetector(rows, cols,alphabetVariance),_forgettingFactor(forgettingFactor),_invForgettingFactor(1.0/forgettingFactor),_nSymbolsToBeDetected(nSymbolsToBeDetected),_alphaPowerSumNow(1.0),_alphaPowerSumPrevious(1.0),_alphaPower(1.0),_g(_channelMatrixRows),_invRtilde(LaGenMatDouble::eye(_channelMatrixRows)),_filter(_channelMatrixRows,_nSymbolsToBeDetected)
// auxiliary
,_identityL(LaGenMatDouble::eye(_channelMatrixRows)),_gObservations(_channelMatrixRows,_channelMatrixRows),_identityMinusgObservations(_channelMatrixRows,_channelMatrixRows),_auxInvRtilde(_channelMatrixRows,_channelMatrixRows),_E(LaGenMatDouble::zeros(_channelMatrixCols,nSymbolsToBeDetected)),_varianceInvRtildeChannelMatrix(_channelMatrixRows,_channelMatrixCols)
{
	tRange rowsRange(_channelMatrixCols-_nSymbolsToBeDetected,_channelMatrixCols-1);
	tRange colsRange(0,_nSymbolsToBeDetected-1);

	_E(rowsRange,colsRange).inject(LaGenMatDouble::eye(_nSymbolsToBeDetected));
}

void RMMSEDetector::StateStep(tVector observations)
{
	if(observations.size()!= _channelMatrixRows)
		throw RuntimeException("observations vector dimensions are wrong.");

	// _g = _invRtilde*observations
	Blas_Mat_Vec_Mult(_invRtilde,observations,_g);

	_alphaPowerSumFactor = _alphaPowerSumNow/(_alphaPowerSumPrevious*_forgettingFactor);

	// _g = _g / (_alphaPowerSumNow + observations'*_invRtilde*observations*_alphaPowerSumFactor
	_g *= 1.0/(_alphaPowerSumNow + Blas_Dot_Prod(observations,_g)*_alphaPowerSumFactor);

	// _gObservations = _g*observations
	Util::Mult(_g,observations,_gObservations);

	// _identityMinusgObservations = _identityL - _alphaPowerSumFactor*_gObservations
	Util::Add(_identityL,_gObservations,_identityMinusgObservations,1.0,-_alphaPowerSumFactor);

	// _auxInvRtilde = _alphaPowerSumFactor*_identityMinusgObservations*_invRtilde
	Blas_Mat_Mat_Mult(_identityMinusgObservations,_invRtilde,_auxInvRtilde,_alphaPowerSumFactor);

	_invRtilde = _auxInvRtilde;

	_alphaPowerSumPrevious = _alphaPowerSumNow;
	_alphaPower *= _forgettingFactor;
	_alphaPowerSumNow = _alphaPowerSumPrevious + _alphaPower;
}

tVector RMMSEDetector::Detect(tVector observations, tMatrix channelMatrix,const tMatrix &noiseCovariance)
{
	if(observations.size()!= _channelMatrixRows || channelMatrix.cols()!=_channelMatrixCols || channelMatrix.rows()!=_channelMatrixRows)
	{
		cout << "channel matrix:" << endl << channelMatrix << endl;
		cout << "channelMatrix.rows(): " << channelMatrix.rows() << " channelMatrix.cols(): " << channelMatrix.cols() << endl;
		throw RuntimeException("RMMSEDetector::Detect: observations vector or channel matrix dimensions are wrong.");
	}

	// the inverse of the observations correlation matrix is updated
	this->StateStep(observations);

	// _varianceInvRtildeChannelMatrix = _alphabetVariance*_invRtilde*channelMatrix
	Blas_Mat_Mat_Mult(_invRtilde,channelMatrix,_varianceInvRtildeChannelMatrix,_alphabetVariance);

	// _filter = _varianceInvRtildeChannelMatrix*_E
	Blas_Mat_Mat_Mult(_varianceInvRtildeChannelMatrix,_E,_filter);

    tVector res(_nSymbolsToBeDetected);

	// _softEstimations = _filter'*observations
	Blas_Mat_Trans_Vec_Mult(_filter,observations,res);

    return res;
}

RMMSEDetector *RMMSEDetector::Clone()
{
	return new RMMSEDetector(*this);
}
