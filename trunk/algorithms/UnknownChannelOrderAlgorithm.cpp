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

UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm(string name, Alphabet alphabet, int L, int N, int frameLength,vector<ChannelMatrixEstimator *> channelEstimators,tMatrix preamble,int iFirstObservation): UnknownChannelAlgorithm(name, alphabet, L, N, frameLength),_channelEstimators(channelEstimators.size()),_candidateOrders( channelEstimators.size()),_maxOrder(-1),_iFirstObservation(iFirstObservation),_preamble(preamble)
{
    for(uint i=0;i<channelEstimators.size();i++)
    {
        // the memory of this estimator is obtained from the number of columns of the channel matrix estimator and the algorithm parameter N
        if((channelEstimators[i]->cols() % _N) !=0)
            throw RuntimeException("UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm: the number of columns of (at least) one of the estimators is not coherent with the number of transmitting antennas (N).");

        // the memory associated with this channel estimator is stored in the corresponding position of _candidateOrders
         _candidateOrders[i] = channelEstimators[i]->cols() / _N;

        // the maximun of the channel matrix estimator orders is obtained
        if(_candidateOrders[i]>_maxOrder)
            _maxOrder = _candidateOrders[i];

        // each of the channel matrix estimators is copied into _channelEstimators
        _channelEstimators[i] = channelEstimators[i]->Clone();
    }

    if(_preamble.cols()<(_maxOrder-1))
        throw RuntimeException("UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm: the number of vectors contained in the preamble are not enough for the maximun channel order of all the channel matrix estimators.");

    // a vector that associate a channel order with its corresponding index is generated
    _channelOrder2index = new int[_maxOrder+1];
    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _channelOrder2index[_candidateOrders[iChannelOrder]] = iChannelOrder;

    _NmaxOrder = _N*_maxOrder;

	_channelOrderAPPs = LaGenMatDouble::zeros(_candidateOrders.size(),_K+_maxOrder-1);
	_channelOrderAPPs(tRange(),tRange(0,_preamble.cols()-1)) = 1.0/double(_candidateOrders.size());

#ifdef DEBUG
    cout << "UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm: _candidateOrders.size " << _candidateOrders.size() << endl;
#endif
}


UnknownChannelOrderAlgorithm::~UnknownChannelOrderAlgorithm()
{
    for(uint i=0;i<_channelEstimators.size();i++)
        delete _channelEstimators[i];

    delete[] _channelOrder2index;
}

vector<vector<tMatrix> > UnknownChannelOrderAlgorithm::EstimateChannelFromTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
    tMatrix sequenceToProcess = Util::append(_preamble,trainingSequence);

    if(observations.cols() < (_iFirstObservation+trainingSequence.cols()))
        throw RuntimeException("UnknownChannelOrderAlgorithm::EstimateChannelFromTrainingSequence: Insufficient number of observations.");

    vector<vector<tMatrix> > estimatedMatrices(_candidateOrders.size());

 	// selects all the rows from a symbols matrix
    tRange rAllSymbolRows(0,_N-1);

    uint iOrder;
    for(int i=_iFirstObservation;i<_iFirstObservation+trainingSequence.cols();i++)
    {
        for(iOrder=0;iOrder<_candidateOrders.size();iOrder++)
        {
            tRange mColumns(_preamble.cols()+i-_iFirstObservation-_candidateOrders[iOrder]+1,_preamble.cols()+i-_iFirstObservation);

            estimatedMatrices[iOrder].push_back(_channelEstimators[iOrder]->nextMatrix(observations.col(i),sequenceToProcess(rAllSymbolRows,mColumns),noiseVariances[i]));
        }
    }

    return estimatedMatrices;
}

