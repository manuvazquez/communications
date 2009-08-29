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
#include "KnownSymbolsCMEapplyingAlgorithm.h"

KnownSymbolsCMEapplyingAlgorithm::KnownSymbolsCMEapplyingAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, MatrixXd preamble, const MatrixXd &symbolVectors): CMEapplyingAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, preamble),_symbolVectors(symbolVectors)
{
}

std::vector<MatrixXd> KnownSymbolsCMEapplyingAlgorithm::estimatedChannelMatricesForChannelOrder(uint iChannelOrder,const MatrixXd &observations,const vector<double> &noiseVariances,const MatrixXd& trainingSequence)
{
    vector<MatrixXd> estimatedChannelMatrices;

    // channel estimation
    for(int iSymbolVector=_preamble.cols();iSymbolVector<_iLastSymbolVectorToBeDetected;iSymbolVector++)
        estimatedChannelMatrices.push_back(_channelEstimators[iChannelOrder]->nextMatrix(observations.col(iSymbolVector),_symbolVectors.block(0,iSymbolVector-_candidateOrders[iChannelOrder]+1,_nInputs,_candidateOrders[iChannelOrder]),noiseVariances[iSymbolVector]));

    return estimatedChannelMatrices;
}

MatrixXd KnownSymbolsCMEapplyingAlgorithm::detectedSymbolsForChannelOrder(uint iChannelOrder,const MatrixXd &observations,const vector<double> &noiseVariances,const MatrixXd& trainingSequence)
{
    return _symbolVectors.block(0,_preamble.cols(),_nInputs,_iLastSymbolVectorToBeDetected-_preamble.cols());
}
