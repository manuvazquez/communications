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

// #define DEBUG

CDMAKalmanEstimator::CDMAKalmanEstimator(tMatrix initialEstimation, int N, const vector<double> &ARprocCoeffs, double ARprocVar): ChannelMatrixEstimator(initialEstimation, N),_stateVectorLength(initialEstimation.cols()*ARprocCoeffs.size())
{
    if(initialEstimation.rows()!=1)
        throw RuntimeException("CDMAKalmanEstimator::CDMAKalmanEstimator: initial channel estimation is not a row vector.");

    if(initialEstimation.cols()!=N)
        throw RuntimeException("CDMAKalmanEstimator::CDMAKalmanEstimator: channel is not flat.");

    // stateTransitionMatrix is a blockwise matrix that represents the state transition matrix
    tMatrix stateTransitionMatrix = LaGenMatDouble::zeros(_stateVectorLength,_stateVectorLength);
    
    uint i,j,iBlockCol,iBlockRow;
    
    // state equation matrix is set
    for(iBlockRow=0;iBlockRow<ARprocCoeffs.size()-1;iBlockRow++)
        for(iBlockCol=iBlockRow+1;iBlockCol<ARprocCoeffs.size();iBlockCol++)
            for(j=iBlockCol*initialEstimation.cols(),i=0;i<initialEstimation.cols();j++,i++)
                stateTransitionMatrix(iBlockRow*initialEstimation.cols()+i,j) = 1.0;
    
    for(iBlockCol=0;iBlockCol<ARprocCoeffs.size();iBlockCol++)
        for(j=iBlockCol*initialEstimation.cols(),i=0;i<initialEstimation.cols();j++,i++)
            stateTransitionMatrix((ARprocCoeffs.size()-1)*initialEstimation.cols()+i,j) = ARprocCoeffs[ARprocCoeffs.size()-1-iBlockCol];
    
    // covariance matrix is set
    tMatrix stateEquationCovariance = LaGenMatDouble::zeros(_stateVectorLength,_stateVectorLength);
    
    for(i=_stateVectorLength-initialEstimation.cols();i<_stateVectorLength;i++)
        stateEquationCovariance(i,i) = ARprocVar;

#ifdef DEBUG
    cout << "AR process coefficients ";
    Util::Print(ARprocCoeffs);
    cout << "stateTransitionMatrix = " << endl << stateTransitionMatrix;
    cout << "stateEquationCovariance = " << endl << stateEquationCovariance;
    getchar();
#endif

    tMatrix initialMean = Util::append(initialEstimation,initialEstimation);
    tMatrix initialCovariance = LaGenMatDouble::eye(_stateVectorLength);

    _kalmanFilter = new KalmanFilter(stateTransitionMatrix,stateEquationCovariance,initialMean,initialCovariance,-1);
}


CDMAKalmanEstimator::CDMAKalmanEstimator(int N): ChannelMatrixEstimator(N)
{
}


CDMAKalmanEstimator::~CDMAKalmanEstimator()
{
    delete _kalmanFilter;
}


ChannelMatrixEstimator* CDMAKalmanEstimator::Clone() const
{
}

double CDMAKalmanEstimator::likelihood(const tVector& observations, const tMatrix symbolsMatrix, double noiseVariance)
{
    return ChannelMatrixEstimator::likelihood(observations, symbolsMatrix, noiseVariance);
}

tMatrix CDMAKalmanEstimator::nextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
}

