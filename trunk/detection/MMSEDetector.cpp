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
#include "MMSEDetector.h"

MMSEDetector::MMSEDetector(int rows, int cols, double alphabetVariance,int nSymbolsToBeDetected): LinearDetector(rows, cols, alphabetVariance),_nSymbolsToBeDetected(nSymbolsToBeDetected),_filter(_channelMatrixRows,_channelMatrixCols)
,_alphabetVarianceChannelMatrixChannelMatrixTrans(rows,rows),_Rx(rows,rows),_piv(rows),_softEstimations(cols),_rNsimbolsDetected(cols-nSymbolsToBeDetected,cols-1),_rAllChannelMatrixRows(0,rows-1)
{
}

MMSEDetector *MMSEDetector::Clone()
{
	return new MMSEDetector(*this);
}

tMatrix MMSEDetector::ComputedFilter()
{
	return _filter(_rAllChannelMatrixRows,_rNsimbolsDetected);
}

tVector MMSEDetector::Detect(tVector observations, tMatrix channelMatrix, const tMatrix& noiseCovariance)
{
	// _alphabetVarianceChannelMatrixChannelMatrixTrans = _alphabetVariance*channelMatrix*channelMatrix^T
	Blas_Mat_Mat_Trans_Mult(channelMatrix,channelMatrix,_alphabetVarianceChannelMatrixChannelMatrixTrans,_alphabetVariance);

	// _Rx = _alphabetVarianceChannelMatrixChannelMatrixTrans + noiseCovariance
	Util::Add(_alphabetVarianceChannelMatrixChannelMatrixTrans,noiseCovariance,_Rx);

	// _invRx = inverse(_Rx)
	LUFactorizeIP(_Rx,_piv);
	LaLUInverseIP(_Rx,_piv);

	// _filter = _Rx*channelMatrix*_alphabetVariance
	Blas_Mat_Mat_Mult(_Rx,channelMatrix,_filter,_alphabetVariance);

	// _softEstimations = _filter'*observations
	Blas_Mat_Trans_Vec_Mult(_filter,observations,_softEstimations);

	return _softEstimations(_rNsimbolsDetected);
}


