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
#include "KnownSymbolsKalmanBasedChannelEstimatorAlgorithm.h"

KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::KnownSymbolsKalmanBasedChannelEstimatorAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m,ChannelMatrixEstimator* channelEstimator, MatrixXd preamble,const MatrixXd &symbolVectors): KnownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected,m, channelEstimator, preamble),_symbolVectors(symbolVectors)
{
}

void KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::run(MatrixXd observations,vector<double> noiseVariances)
{
    _estimatedChannelMatrices.reserve(_iLastSymbolVectorToBeDetected-_preamble.cols());

    for(int iSymbolVector=_preamble.cols();iSymbolVector<_iLastSymbolVectorToBeDetected;iSymbolVector++)
        _estimatedChannelMatrices.push_back(_channelEstimator->nextMatrix(observations.col(iSymbolVector),_symbolVectors.block(0,iSymbolVector-_channelOrder+1,_nInputs,_channelOrder),noiseVariances[iSymbolVector]));
}

void KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::run(MatrixXd observations,vector<double> noiseVariances,MatrixXd trainingSequence)
{
    run(observations,noiseVariances);
}

MatrixXd KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::getDetectedSymbolVectors()
{
    MatrixXd aux(1,1);
    aux.resize(0,0);
    return aux;
}

vector<MatrixXd> KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::getEstimatedChannelMatrices()
{
    return _estimatedChannelMatrices;
}
