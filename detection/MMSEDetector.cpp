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

MMSEDetector::MMSEDetector(int rows, int cols, double alphabetVariance,uint nSymbolsToBeDetected): LinearDetector(rows, cols, alphabetVariance),_nSymbolsToBeDetected(nSymbolsToBeDetected),_detectionStart(cols-nSymbolsToBeDetected)
{
}

MMSEDetector::MMSEDetector(int rows, int cols, double alphabetVariance,uint nSymbolsToBeDetected,uint startingFrom): LinearDetector(rows, cols, alphabetVariance),_nSymbolsToBeDetected(nSymbolsToBeDetected),_detectionStart(startingFrom)
{
	if(_detectionStart+_nSymbolsToBeDetected>_channelMatrixCols)
		throw RuntimeException("MMSEDetector::MMSEDetector: nSymbolsToBeDetected, startingFrom or both parameters are wrong.");
}

MMSEDetector *MMSEDetector::clone()
{
	return new MMSEDetector(*this);
}

MatrixXd MMSEDetector::computedFilter_eigen()
{
    return _filter_eigen.block(0,_detectionStart,_channelMatrixRows,_nSymbolsToBeDetected);
}

VectorXd MMSEDetector::detect(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance)
{
    MatrixXd _Rx_eigen = noiseCovariance + _alphabetVariance*channelMatrix*channelMatrix.transpose();

    _filter_eigen = _Rx_eigen.inverse()*channelMatrix*_alphabetVariance;

    VectorXd softEstimations = _filter_eigen.transpose()*observations;

    // required for nthSymbolVariance computing
    _channelMatrix_eigen = channelMatrix;

    return softEstimations.segment(_detectionStart,_nSymbolsToBeDetected);
}

double MMSEDetector::nthSymbolVariance(int n,double noiseVariance)
{
    return nthSymbolGain(n)*(1.0-nthSymbolGain(n));
}

double MMSEDetector::nthSymbolGain(int n) const
{
    return _filter_eigen.col(_detectionStart+n).dot(_channelMatrix_eigen.col(_detectionStart+n));   
}
