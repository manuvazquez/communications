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

KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::KnownSymbolsKalmanBasedChannelEstimatorAlgorithm(string name, Alphabet alphabet,int L,int N, int frameLength,int m,ChannelMatrixEstimator* channelEstimator, tMatrix preamble,const tMatrix &symbolVectors): KnownChannelOrderAlgorithm(name, alphabet, L, N, frameLength,m, channelEstimator, preamble),_symbolVectors(symbolVectors)
{
}

KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::~KnownSymbolsKalmanBasedChannelEstimatorAlgorithm()
{
}

void KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    _estimatedChannelMatrices.reserve(_K-_preamble.cols());

    tRange rAllSymbolRows(0,_nInputs-1);

    for(int iSymbolVector=_preamble.cols();iSymbolVector<_K;iSymbolVector++)
    {
        _estimatedChannelMatrices.push_back( _channelEstimator->nextMatrix(observations.col(iSymbolVector),_symbolVectors(rAllSymbolRows,tRange(iSymbolVector-_channelOrder+1,iSymbolVector)),noiseVariances[iSymbolVector]));
    }
}

void KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::Run(tMatrix observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
    Run(observations,noiseVariances);
}

tMatrix KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::getDetectedSymbolVectors()
{
    return tMatrix(0,0);
}

vector<tMatrix> KnownSymbolsKalmanBasedChannelEstimatorAlgorithm::GetEstimatedChannelMatrices()
{
    return _estimatedChannelMatrices;
}
