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
#include "UnknownChannelOrderAlgorithm.h"

// #define DEBUG

UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected,vector<ChannelMatrixEstimator *> channelEstimators,MatrixXd preamble,int iFirstObservation): UnknownChannelAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected),_channelEstimators(channelEstimators.size()),_candidateOrders( channelEstimators.size()),_maxOrder(-1),_iFirstObservation(iFirstObservation),_preamble(preamble)
{
    for(uint i=0;i<channelEstimators.size();i++)
    {
        // the memory of this estimator is obtained from the number of columns of the channel matrix estimator and the algorithm parameter N
        if((channelEstimators[i]->cols() % _nInputs) !=0)
            throw RuntimeException("UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm: the number of columns of (at least) one of the estimators is not coherent with the number of transmitting antennas (N).");

        // the memory associated with this channel estimator is stored in the corresponding position of _candidateOrders
         _candidateOrders[i] = channelEstimators[i]->cols() / _nInputs;

        // the maximun of the channel matrix estimator orders is obtained
        if(_candidateOrders[i]>_maxOrder)
            _maxOrder = _candidateOrders[i];

        // each of the channel matrix estimators is copied into _channelEstimators
        _channelEstimators[i] = channelEstimators[i]->clone();
    }

    if(_preamble.cols()<(_maxOrder-1))
        throw RuntimeException("UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm: the number of vectors contained in the preamble are not enough for the maximun channel order of all the channel matrix estimators.");

    // a vector that associate a channel order with its corresponding index is generated
    _channelOrder2index = new int[_maxOrder+1];
    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _channelOrder2index[_candidateOrders[iChannelOrder]] = iChannelOrder;

    _nInputsXmaxChannelOrder = _nInputs*_maxOrder;

    _channelOrderAPPs = MatrixXd::Zero(_candidateOrders.size(),_iLastSymbolVectorToBeDetected+_maxOrder-1);
    _channelOrderAPPs.block(0,0,_candidateOrders.size(),_preamble.cols()).setConstant(1.0/double(_candidateOrders.size()));
}


UnknownChannelOrderAlgorithm::~UnknownChannelOrderAlgorithm()
{
    for(uint i=0;i<_channelEstimators.size();i++)
        delete _channelEstimators[i];

    delete[] _channelOrder2index;
}
