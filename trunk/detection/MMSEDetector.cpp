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

// #define DEBUG

MMSEDetector::MMSEDetector(int rows, int cols, double alphabetVariance,int nSymbolsToBeDetected): LinearDetector(rows, cols, alphabetVariance),_nSymbolsToBeDetected(nSymbolsToBeDetected),_detectionStart(cols-nSymbolsToBeDetected),_filter(_channelMatrixRows,_channelMatrixCols)
{
}

MMSEDetector::MMSEDetector(int rows, int cols, double alphabetVariance,int nSymbolsToBeDetected,int startingFrom): LinearDetector(rows, cols, alphabetVariance),_nSymbolsToBeDetected(nSymbolsToBeDetected),_detectionStart(startingFrom),_filter(_channelMatrixRows,_channelMatrixCols)
{
	if(_detectionStart+_nSymbolsToBeDetected>_channelMatrixCols)
		throw RuntimeException("MMSEDetector::MMSEDetector: nSymbolsToBeDetected, startingFrom or both parameters are wrong.");
}

MMSEDetector *MMSEDetector::Clone()
{
	return new MMSEDetector(*this);
}

tMatrix MMSEDetector::ComputedFilter()
{
	tRange rNsimbolsDetected(_detectionStart,_detectionStart+_nSymbolsToBeDetected-1);
	return _filter(tRange(),rNsimbolsDetected);
// 	return _filter;
}

tVector MMSEDetector::Detect(tVector observations, tMatrix channelMatrix, const tMatrix& noiseCovariance)
{
	tMatrix alphabetVarianceChannelMatrixChannelMatrixTrans(_channelMatrixRows,_channelMatrixRows);
	// alphabetVarianceChannelMatrixChannelMatrixTrans = _alphabetVariance*channelMatrix*channelMatrix^T
	Blas_Mat_Mat_Trans_Mult(channelMatrix,channelMatrix,alphabetVarianceChannelMatrixChannelMatrixTrans,_alphabetVariance);

	tMatrix Rx(_channelMatrixRows,_channelMatrixRows);
	// Rx = alphabetVarianceChannelMatrixChannelMatrixTrans + noiseCovariance
	Util::add(alphabetVarianceChannelMatrixChannelMatrixTrans,noiseCovariance,Rx);

	_Rx = Rx;

	tLongIntVector piv(_channelMatrixRows);
	// _invRx = inverse(Rx)
	LUFactorizeIP(Rx,piv);
	LaLUInverseIP(Rx,piv);

	// _filter = Rx*channelMatrix*_alphabetVariance
	Blas_Mat_Mat_Mult(Rx,channelMatrix,_filter,_alphabetVariance);

	tVector softEstimations(_channelMatrixCols);
	// _softEstimations = _filter'*observations
	Blas_Mat_Trans_Vec_Mult(_filter,observations,softEstimations);

	// ----------------- required for nthSymbolVariance computing -------------------
	_channelMatrix = channelMatrix;

	// ------------------------------------------------------------------------------

	// we only keep the columns and soft estimations we are interested in
	tRange rNsimbolsDetected(_detectionStart,_detectionStart+_nSymbolsToBeDetected-1);

// 	_filter = _filter(tRange(),rNsimbolsDetected);


#ifdef DEBUG
	tMatrix aux(_channelMatrixCols,_channelMatrixCols);
	Blas_Mat_Trans_Mat_Mult(_filter,channelMatrix,aux);
	cout << "filtro * canal" << endl << aux;
#endif

// 	return softEstimations(_rNsimbolsDetected);

	return softEstimations(rNsimbolsDetected);
}

double MMSEDetector::nthSymbolVariance(int n)
{
#ifdef DEBUG
	cout << "mu = " << Blas_Dot_Prod(_filter.col(_detectionStart+n),_channelMatrix.col(_detectionStart+n)) << endl;
#endif
	tVector Rxf(_channelMatrixRows);

#ifdef DEBUG
	cout << "Lo que devolvía antes " << (1.0 - Blas_Dot_Prod(_filter.col(_detectionStart+n),_channelMatrix.col(_detectionStart+n))) << endl;
#endif

	// Rxf = _Rx * filter.col(_channelMatrixCols-_nSymbolsToBeDetected+n)
	Blas_Mat_Vec_Mult(_Rx,_filter.col(_detectionStart+n),Rxf);

#ifdef DEBUG
	cout << "Lo que voy a devolver ahora: " << Blas_Dot_Prod(_filter.col(_detectionStart+n),Rxf) - pow(Blas_Dot_Prod(_filter.col(_detectionStart+n),_channelMatrix.col(_detectionStart+n)),2.0) << endl;
#endif

// 	return Blas_Dot_Prod(_filter.col(_detectionStart+n),Rxf) - pow(Blas_Dot_Prod(_filter.col(_detectionStart+n),_channelMatrix.col(_detectionStart+n)),2.0);
	return (1.0 - Blas_Dot_Prod(_filter.col(_channelMatrixCols-_nSymbolsToBeDetected+n),_channelMatrix.col(_channelMatrixCols-_nSymbolsToBeDetected+n)));
}

// double MMSEDetector::nthSymbolGain(int n) const
// {
// 	return Blas_Dot_Prod(_filter.col(_detectionStart+n),_channelMatrix.col(_detectionStart+n));
// }
