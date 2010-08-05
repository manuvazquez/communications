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

// #define PRINT_INFO

CDMAKalmanEstimator::CDMAKalmanEstimator(const MatrixXd& initialEstimation, const MatrixXd& variances, vector< double > ARcoefficients, double ARvariance, const MatrixXd &spreadingCodes): KalmanEstimator(initialEstimation, variances, spreadingCodes.cols(), ARcoefficients, ARvariance),_spreadingCodes(spreadingCodes)
{
    if(spreadingCodes.cols()!=_nInputs)
        throw RuntimeException("CDMAKalmanEstimator::CDMAKalmanEstimator: the number of spreading codes doesn't match the number of users.");
    
    _nOutputs = _spreadingCodes.rows();
}


CDMAKalmanEstimator::CDMAKalmanEstimator(const CDMAKalmanEstimator& cdmaKalmanEstimator):KalmanEstimator(cdmaKalmanEstimator),_spreadingCodes(cdmaKalmanEstimator._spreadingCodes)
{
}

CDMAKalmanEstimator* CDMAKalmanEstimator::clone() const
{
    return new CDMAKalmanEstimator(*this);
}

MatrixXd CDMAKalmanEstimator::buildMeasurementMatrix(const VectorXd& symbolsVector)
{
    if(symbolsVector.size()!=_nInputs)
        throw RuntimeException("CDMAKalmanEstimator::buildMeasurementMatrix: symbols vector length is wrong.");
            
    MatrixXd CS = MatrixXd::Zero(_nOutputs,_nInputs);
    
    for(int i=0;i<_nOutputs;i++)
        for(int j=0;j<_nInputs;j++)
            CS(i,j) = _spreadingCodes(i,j)*symbolsVector(j);            
            
#ifdef PRINT_INFO
    cout << "CDMAKalmanEstimator::BuildFfromSymbolsMatrix: measurement matrix built" << endl << CS;
#endif            
            
    return CS;
}

MatrixXd CDMAKalmanEstimator::sampleFromPredictive() const
{
    MatrixXd sampledChannelMatrix = KalmanEstimator::sampleFromPredictive();
    
    if(sampledChannelMatrix.rows()!=1)
        throw RuntimeException("CDMAKalmanEstimator::sampleFromPredictive_eigen: sampled channel matrix is not a row vector.");
    
    MatrixXd spreadingCodesXsampledChannelMatrix = _spreadingCodes;
    for(int i=0;i<_nOutputs;i++)
        for(int j=0;j<_nInputs;j++)
            spreadingCodesXsampledChannelMatrix(i,j) *= sampledChannelMatrix(0,j);
            
    return spreadingCodesXsampledChannelMatrix;
}

MatrixXd CDMAKalmanEstimator::lastEstimatedChannelMatrix() const
{
    MatrixXd sampledChannelMatrix = KalmanEstimator::lastEstimatedChannelMatrix();
    
    if(sampledChannelMatrix.rows()!=1)
        throw RuntimeException("CDMAKalmanEstimator::lastEstimatedChannelMatrix_eigen: sampled channel matrix is not a row vector.");
    
    MatrixXd spreadingCodesXsampledChannelMatrix = _spreadingCodes;
    for(int i=0;i<_nOutputs;i++)
        for(int j=0;j<_nInputs;j++)
            spreadingCodesXsampledChannelMatrix(i,j) *= sampledChannelMatrix(0,j);
            
    return spreadingCodesXsampledChannelMatrix;
}

MatrixXd CDMAKalmanEstimator::getPredictive() const
{
    MatrixXd sampledChannelMatrix = KalmanEstimator::getPredictive();
    
    if(sampledChannelMatrix.rows()!=1)
        throw RuntimeException("CDMAKalmanEstimator::getPredictive: sampled channel matrix is not a row vector.");
    
    MatrixXd spreadingCodesXsampledChannelMatrix = _spreadingCodes;
    for(int i=0;i<_nOutputs;i++)
        for(int j=0;j<_nInputs;j++)
            spreadingCodesXsampledChannelMatrix(i,j) *= sampledChannelMatrix(0,j);
            
    return spreadingCodesXsampledChannelMatrix;
}