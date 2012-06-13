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

RMMSEDetector::RMMSEDetector(uint rows, uint cols,double alphabetVariance,double forgettingFactor,uint nSymbolsToBeDetected): LinearDetector(rows, cols,alphabetVariance),_forgettingFactor(forgettingFactor),_invForgettingFactor(1.0/forgettingFactor),_nSymbolsToBeDetected(nSymbolsToBeDetected),_alphaPowerSumNow(1.0),_alphaPowerSumPrevious(1.0),_alphaPower(1.0),_invRtilde(MatrixXd::Identity(rows,rows)),_E(MatrixXd::Zero(cols,nSymbolsToBeDetected))
{
    _E.block(_channelMatrixCols-_nSymbolsToBeDetected,0,_nSymbolsToBeDetected,_nSymbolsToBeDetected) = MatrixXd::Identity(_nSymbolsToBeDetected,_nSymbolsToBeDetected);
}

void RMMSEDetector::stateStep(const VectorXd &observations)
{
    if(observations.size()!= _channelMatrixRows)
    {
        cout << "observations.size() = " << observations.size() << " _channelMatrixRows = " << _channelMatrixRows << endl;
        throw RuntimeException("RMMSEDetector::StateStep: observations vector dimensions are wrong.");
    }

    _alphaPowerSumFactor = _alphaPowerSumNow/(_alphaPowerSumPrevious*_forgettingFactor);

    _g = _invRtilde*observations/(_alphaPowerSumNow + observations.dot(_invRtilde*observations)*_alphaPowerSumFactor);
    
    _invRtilde = _alphaPowerSumFactor*(MatrixXd::Identity(_channelMatrixRows,_channelMatrixRows) - _alphaPowerSumFactor*_g*observations.transpose())*_invRtilde;    

    _alphaPowerSumPrevious = _alphaPowerSumNow;
    _alphaPower *= _forgettingFactor;
    _alphaPowerSumNow = _alphaPowerSumPrevious + _alphaPower;
}

VectorXd RMMSEDetector::detect(const VectorXd &observations, const MatrixXd &channelMatrix,const MatrixXd &noiseCovariance)
{
    if(observations.size()!= _channelMatrixRows || channelMatrix.cols()!=_channelMatrixCols || channelMatrix.rows()!=_channelMatrixRows)
    {
        cout << "channel matrix:" << endl << channelMatrix << endl;
        cout << "channelMatrix.rows(): " << channelMatrix.rows() << " channelMatrix.cols(): " << channelMatrix.cols() << endl;
        throw RuntimeException("RMMSEDetector::Detect: observations vector or channel matrix dimensions are wrong.");
    }

    // the inverse of the observations correlation matrix is updated
    stateStep(observations);

    _filter = _alphabetVariance*_invRtilde*channelMatrix*_E;

    // required for nthSymbolVariance computing
    _alphabetVarianceChannelMatrixChannelMatrixTransPlusNoiseCovariance = _alphabetVariance*channelMatrix*channelMatrix.transpose()+noiseCovariance;
    _channelMatrix = channelMatrix;

    return _filter.transpose()*observations;
}

RMMSEDetector *RMMSEDetector::clone()
{
	return new RMMSEDetector(*this);
}

double RMMSEDetector::nthSymbolVariance(uint n,double noiseVariance) const
{
    return _alphabetVariance*(1.0 - 2.0*_filter.col(n).dot(_channelMatrix.col(_channelMatrixCols-_nSymbolsToBeDetected+n))) + _filter.col(n).dot(_alphabetVarianceChannelMatrixChannelMatrixTransPlusNoiseCovariance*_filter.col(n));
}
