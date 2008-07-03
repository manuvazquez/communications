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
#include "KnownChannelOrderAlgorithm.h"

KnownChannelOrderAlgorithm::KnownChannelOrderAlgorithm(string name, Alphabet alphabet,int L,int N, int K,int m, ChannelMatrixEstimator *channelEstimator,tMatrix preamble): UnknownChannelAlgorithm(name, alphabet,L,N,K),_channelEstimator(channelEstimator->Clone()),_m(m),_Nm(channelEstimator->Cols()),_preamble(preamble)
{
    if(_m!=(_Nm/_N))
        throw RuntimeException("KnownChannelOrderAlgorithm::KnownChannelOrderAlgorithm: the channel order parameter is not coherent with the channel estimator.");
}

KnownChannelOrderAlgorithm::KnownChannelOrderAlgorithm(string name, Alphabet alphabet,int L,int N, int K,int m,tMatrix preamble): UnknownChannelAlgorithm(name, alphabet,L,N,K),_channelEstimator(NULL),_m(m),_Nm(_N*m),_preamble(preamble)
{
}

KnownChannelOrderAlgorithm::~ KnownChannelOrderAlgorithm()
{
	delete _channelEstimator;
}

// vector<tMatrix> KnownChannelOrderAlgorithm::EstimateChannelFromTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence,ChannelMatrixEstimator *channelEstimator)
// {
// 	tMatrix toProcessSequence = Util::Append(_preamble,trainingSequence);
// 	int lengthToProcessSequence = toProcessSequence.cols();
//
// 	if(observations.cols()<lengthToProcessSequence)
// 		throw RuntimeException("KnownChannelOrderAlgorithm::EstimateChannelFromTrainingSequence: insufficient number of observations.");
//
// 	vector<tMatrix> estimatedMatrices(trainingSequence.cols());
//
// 	// selects all the rows from a symbols matrix
// 	tRange allSymbolRows(0,_N-1);
//
// 	for(int i=_preamble.cols();i<lengthToProcessSequence;i++)
// 	{
// 		tRange mColumns(i-_m+1,i);
// 		estimatedMatrices[i-_preamble.cols()] = channelEstimator->NextMatrix(observations.col(i),toProcessSequence(allSymbolRows,mColumns),noiseVariances[i]);
// 	}
// 	return estimatedMatrices;
// }
