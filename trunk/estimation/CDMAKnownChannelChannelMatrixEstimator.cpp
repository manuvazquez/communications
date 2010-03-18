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
#include "CDMAKnownChannelChannelMatrixEstimator.h"

// #define DEBUG

CDMAKnownChannelChannelMatrixEstimator::CDMAKnownChannelChannelMatrixEstimator(const MIMOChannel *channel, int iFirstChannelMatrix, int N, const MatrixXd &spreadingCodes): KnownChannelChannelMatrixEstimator(channel, iFirstChannelMatrix, N),_spreadingCodes(spreadingCodes)
{
    if(_channel->at(iFirstChannelMatrix).rows()!=1)
        throw RuntimeException("CDMAKnownChannelChannelMatrixEstimator::CDMAKnownChannelChannelMatrixEstimator: channel matrices don't have a single row.");

    _nOutputs = _spreadingCodes.rows();
}

// eigen
double CDMAKnownChannelChannelMatrixEstimator::likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance)
{
    if(symbolsMatrix.cols()!=1)
        throw RuntimeException("CDMAKnownChannelChannelMatrixEstimator::likelihood: the symbols matrix received should be a column vector.");
    
    MatrixXd channelCoefficientsXsymbols = symbolsMatrix;
        
    for(int i=0;i<_nInputs;i++)
        channelCoefficientsXsymbols(i,0) *= _lastEstimatedChannelMatrix(0,i);
         
    return StatUtil::normalPdf(observations,_spreadingCodes*channelCoefficientsXsymbols,noiseVariance);
}

CDMAKnownChannelChannelMatrixEstimator *CDMAKnownChannelChannelMatrixEstimator::clone() const
{
    return new CDMAKnownChannelChannelMatrixEstimator(*this);
}
