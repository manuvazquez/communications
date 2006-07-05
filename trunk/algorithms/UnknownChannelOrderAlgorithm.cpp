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

UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm(string name, Alphabet alphabet, int L, int N, int K,vector<ChannelMatrixEstimator *> channelEstimators,tMatrix preamble,int firstObservationIndex): UnknownChannelAlgorithm(name, alphabet, L, N, K),_channelEstimators(channelEstimators.size()),_preamble(preamble),_candidateOrders(new int[channelEstimators.size()]),_maxOrder(-1),_firstObservationIndex(firstObservationIndex),_nCandidateOrders(channelEstimators.size())
{
	for(int i=0;i<_nCandidateOrders;i++)
	{
		// the memory of this estimator is obtained from the number of columns of the channel matrix estimator and the algorithm parameter N
		if((channelEstimators[i]->Cols() % _N) !=0)
			throw RuntimeException("UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm: the number of columns of (at least) one of the estimators is not coherent with the number of transmitting antennas (N).");

		// the memory associated with this channel estimator is stored in the corresponding position of _candidateOrders
		 _candidateOrders[i] = channelEstimators[i]->Cols() / _N;

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
	for(int iChannelOrder=0;iChannelOrder<_nCandidateOrders;iChannelOrder++)
		_channelOrder2index[_candidateOrders[iChannelOrder]] = iChannelOrder;
}


UnknownChannelOrderAlgorithm::~UnknownChannelOrderAlgorithm()
{
	for(int i=0;i<_channelEstimators.size();i++)
		delete _channelEstimators[i];

	delete[] _candidateOrders;
	delete[] _channelOrder2index;
}

vector<vector<tMatrix> > UnknownChannelOrderAlgorithm::ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
	tMatrix sequenceToProcess = Util::Append(_preamble,trainingSequence);
	int lengthSequenceToProcess = sequenceToProcess.cols();
	
	if(observations.cols() < (_firstObservationIndex+trainingSequence.cols()))
		throw RuntimeException("UnknownChannelOrderAlgorithm::ProcessTrainingSequence: Insufficient number of observations.");
	
	vector<vector<tMatrix> > estimatedMatrices(_nCandidateOrders);

// 	// selects all the rows from a symbols matrix
	tRange rAllSymbolRows(0,_N-1);

	int iOrder;
	for(int i=_firstObservationIndex;i<_firstObservationIndex+trainingSequence.cols();i++)
	{
		for(iOrder=0;iOrder<_nCandidateOrders;iOrder++)
		{
			tRange mColumns(_preamble.cols()+i-_firstObservationIndex-_candidateOrders[iOrder]+1,_preamble.cols()+i-_firstObservationIndex);
//
			estimatedMatrices[iOrder].push_back(_channelEstimators[iOrder]->NextMatrix(observations.col(i),sequenceToProcess(rAllSymbolRows,mColumns),noiseVariances[i]));
		}
	}

// 	for(int j=0;j<_nCandidateOrders;j++)
// 		cout << _channelEstimators[j]->LastEstimatedChannelMatrix();

	return estimatedMatrices;
}

