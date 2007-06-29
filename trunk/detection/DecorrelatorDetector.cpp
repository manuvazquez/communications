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
#include "DecorrelatorDetector.h"

DecorrelatorDetector::DecorrelatorDetector(int rows, int cols, double alphabetVariance): LinearDetector(rows, cols, alphabetVariance),_filter(_channelMatrixCols,_channelMatrixRows)
{
}

double DecorrelatorDetector::NthSymbolVariance(int n)
{
	return Blas_Dot_Prod(_filter.row(n),_filter.row(n));
}

LinearDetector* DecorrelatorDetector::Clone()
{
	return new DecorrelatorDetector(*this);
}

tVector DecorrelatorDetector::Detect(tVector observations, tMatrix channelMatrix, const tMatrix& noiseCovariance)
{
// 	std::cout << "hola" << std::endl;
	tMatrix channelMatrixChannelMatrixTrans(_channelMatrixCols,_channelMatrixCols);

	// channelMatrixChannelMatrixTrans = _alphabetVariance*channelMatrix*channelMatrix^T
	Blas_Mat_Trans_Mat_Mult(channelMatrix,channelMatrix,channelMatrixChannelMatrixTrans);

	tLongIntVector piv(_channelMatrixRows);

	// channelMatrixChannelMatrixTrans = inverse(channelMatrixChannelMatrixTrans)
	LUFactorizeIP(channelMatrixChannelMatrixTrans,piv);
	LaLUInverseIP(channelMatrixChannelMatrixTrans,piv);

	// _filter = channelMatrixChannelMatrixTrans*channelMatrix*channelMatrix^T
	Blas_Mat_Mat_Trans_Mult(channelMatrixChannelMatrixTrans,channelMatrix,_filter);

	tVector softEstimations(_channelMatrixCols);

	// softEstimations = _filter*observations
	Blas_Mat_Vec_Mult(_filter,observations,softEstimations);

	return softEstimations;
}

