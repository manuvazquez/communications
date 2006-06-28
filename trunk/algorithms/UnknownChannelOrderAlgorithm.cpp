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

UnknownChannelOrderAlgorithm::UnknownChannelOrderAlgorithm(string name, Alphabet alphabet, int L, int N, int K,vector<ChannelMatrixEstimator *> channelEstimators,tMatrix preamble): UnknownChannelAlgorithm(name, alphabet, L, N, K),_channelEstimators(channelEstimators.size()),_preamble(preamble),_candidateOrders(new int[channelEstimators.size()]),_maxOrder(-1)
{
	for(int i=0;i<channelEstimators.size();i++)
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
}


UnknownChannelOrderAlgorithm::~UnknownChannelOrderAlgorithm()
{
	for(int i=0;i<_channelEstimators.size();i++)
		delete _channelEstimators[i];

	delete[] _candidateOrders;
}

vector<vector<tMatrix> > UnknownChannelOrderAlgorithm::ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
// // 	int lengthSequenceToProcess = trainingSequence.cols() + _preamble.cols();
// 	tMatrix toProcessSequence = Util::Append(_preamble,trainingSequence);
// 	int lengthToProcessSequence = toProcessSequence.cols();
// 	
// 	if(observations.cols()<lengthToProcessSequence)
// 		throw RuntimeException("Insufficient number of observations.");
// 	
// 	vector<tMatrix> estimatedMatrices(lengthToProcessSequence);
// 
// 	// selects all the rows from a symbols matrix
// 	tRange allSymbolRows(0,_N-1);
// 
// 	for(int i=_m-1;i<lengthToProcessSequence;i++)
// 	{
// 		tRange mColumns(i-_m+1,i);
// 		estimatedMatrices[i] = _channelEstimator->NextMatrix(observations.col(i),toProcessSequence(allSymbolRows,mColumns),noiseVariances[i]);
// 	}
// 	return estimatedMatrices;
}

