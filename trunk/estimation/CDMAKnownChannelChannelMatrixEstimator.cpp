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

CDMAKnownChannelChannelMatrixEstimator::CDMAKnownChannelChannelMatrixEstimator(const MIMOChannel *channel, int iFirstChannelMatrix, int N, const tMatrix &spreadingCodes): KnownChannelChannelMatrixEstimator(channel, iFirstChannelMatrix, N),_spreadingCodes(spreadingCodes)
{
    if((*_channel)[iFirstChannelMatrix].rows()!=1)
        throw RuntimeException("CDMAKnownChannelChannelMatrixEstimator::CDMAKnownChannelChannelMatrixEstimator: channel matrices don't have a single row.");

    _nOutputs = _spreadingCodes.rows();
}

double CDMAKnownChannelChannelMatrixEstimator::likelihood(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance)
{
    if(symbolsMatrix.cols()!=1)
        throw RuntimeException("CDMAKnownChannelChannelMatrixEstimator::likelihood: the symbols matrix received should be a column vector.");
    
    tMatrix channelCoefficientsXsymbols(symbolsMatrix);
        
    for(int i=0;i<_nInputs;i++)
        channelCoefficientsXsymbols(i,0) *= _lastEstimatedChannelMatrix(0,i);

    tVector withoutNoiseObservations(_nOutputs);    
    
    Blas_Mat_Vec_Mult(_spreadingCodes,channelCoefficientsXsymbols,withoutNoiseObservations);
    
#ifdef DEBUG
//     cout << "observations = " << endl << observations;
    cout << "_lastEstimatedChannelMatrix = " << endl << _lastEstimatedChannelMatrix;
    cout << "symbolsMatrix = " << endl << symbolsMatrix;
//     cout << "channelCoefficientsXsymbols = " << endl << channelCoefficientsXsymbols;
//     cout << "_spreadingCodes = " << endl << _spreadingCodes;
    cout << "withoutNoiseObservations = " << endl << withoutNoiseObservations << endl;
    cout << "likelihood = " << StatUtil::NormalPdf(observations,withoutNoiseObservations,noiseVariance) << endl;
//     getchar();
#endif        
    
    return StatUtil::NormalPdf(observations,withoutNoiseObservations,noiseVariance);
}

CDMAKnownChannelChannelMatrixEstimator *CDMAKnownChannelChannelMatrixEstimator::clone() const
{
    return new CDMAKnownChannelChannelMatrixEstimator(*this);
}
