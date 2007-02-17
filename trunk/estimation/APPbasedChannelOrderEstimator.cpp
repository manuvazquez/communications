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
#include "APPbasedChannelOrderEstimator.h"

APPbasedChannelOrderEstimator::APPbasedChannelOrderEstimator(const tMatrix& preamble,vector<ChannelMatrixEstimator *> channelEstimators): ChannelOrderEstimator(preamble),_channelEstimators(channelEstimators.size())
{
	for(int i=0;i<channelEstimators.size();i++)
		_channelEstimators[i] = channelEstimators[i]->Clone();
}


APPbasedChannelOrderEstimator::~APPbasedChannelOrderEstimator()
{
	for(int i=0;i<_channelEstimators.size();i++)
		delete _channelEstimators[i];
}


vector< double > APPbasedChannelOrderEstimator::ComputeProbabilities(const tMatrix& observations, vector< double > noiseVariances, tMatrix trainingSequence)
{
//     tMatrix sequenceToProcess = Util::Append(_preamble,trainingSequence);
//     int lengthSequenceToProcess = sequenceToProcess.cols();
//
//     if(observations.cols() < (_iFirstObservation+trainingSequence.cols()))
//         throw RuntimeException("UnknownChannelOrderAlgorithm::ProcessTrainingSequence: Insufficient number of observations.");
//
//     vector<vector<tMatrix> > estimatedMatrices(_candidateOrders.size());
//
//  	// selects all the rows from a symbols matrix
//     tRange rAllSymbolRows(0,_N-1);
//
//     int iOrder;
//     for(int i=_iFirstObservation;i<_iFirstObservation+trainingSequence.cols();i++)
//     {
//         for(iOrder=0;iOrder<_candidateOrders.size();iOrder++)
//         {
//             tRange mColumns(_preamble.cols()+i-_iFirstObservation-_candidateOrders[iOrder]+1,_preamble.cols()+i-_iFirstObservation);
// //
//             estimatedMatrices[iOrder].push_back(_channelEstimators[iOrder]->NextMatrix(observations.col(i),sequenceToProcess(rAllSymbolRows,mColumns),noiseVariances[i]));
//         }
//     }
}

