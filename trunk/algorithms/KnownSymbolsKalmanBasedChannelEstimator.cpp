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
#include "KnownSymbolsKalmanBasedChannelEstimator.h"

KnownSymbolsKalmanBasedChannelEstimator::KnownSymbolsKalmanBasedChannelEstimator(string name, Alphabet alphabet, int K, KalmanEstimator* channelEstimator, tMatrix preamble,const tMatrix &symbolVectors): KnownChannelOrderAlgorithm(name, alphabet, K, channelEstimator, preamble),_symbolVectors(symbolVectors)
{
//     cout << "Constructor: " << endl << _symbolVectors << endl; 
//     cout << "Es submatrixview " << _symbolVectors.is_submatrixview() << endl;
}


KnownSymbolsKalmanBasedChannelEstimator::~KnownSymbolsKalmanBasedChannelEstimator()
{
}

void KnownSymbolsKalmanBasedChannelEstimator::Run(tMatrix observations,vector<double> noiseVariances)
{
//         cout << "En Run: " << endl << _symbolVectors << endl;
//     cout << "Es submatrixview " << _symbolVectors.is_submatrixview() << endl;

//     int nSymbolVectors = _symbolVectors.cols();
    _estimatedChannelMatrices.reserve(_K-_m+1);
//     _estimatedChannelMatrices.reserve(nSymbolVectors-_m+1);

    tRange rAllSymbolRows(0,_N-1);

//     for(int iSymbolVector=_m-1;iSymbolVector<nSymbolVectors;iSymbolVector++)
    for(int iSymbolVector=_m-1;iSymbolVector<_K;iSymbolVector++)
    {
        _estimatedChannelMatrices.push_back( _channelEstimator->NextMatrix(observations.col(iSymbolVector),_symbolVectors(rAllSymbolRows,tRange(iSymbolVector-_m+1,iSymbolVector)),noiseVariances[iSymbolVector]));
    }
}

void KnownSymbolsKalmanBasedChannelEstimator::Run(tMatrix observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
    Run(observations,noiseVariances);
}

tMatrix KnownSymbolsKalmanBasedChannelEstimator::GetDetectedSymbolVectors()
{
    return tMatrix(0,0);
}

vector<tMatrix> KnownSymbolsKalmanBasedChannelEstimator::GetEstimatedChannelMatrices()
{
    return _estimatedChannelMatrices;
}
