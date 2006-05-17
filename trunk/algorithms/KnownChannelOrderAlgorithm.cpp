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

KnownChannelOrderAlgorithm::KnownChannelOrderAlgorithm(string name, Alphabet alphabet, ChannelMatrixEstimator *channelEstimator,tMatrix preamble): UnknownChannelAlgorithm(name, alphabet, channelEstimator),_preamble(preamble),_L(channelEstimator->Rows()),_Nm(channelEstimator->Cols())
{
	// if there is no preamble
	if(_preamble.rows()==0)
		_m = 1;
	else
		_m = _preamble.cols() + 1;

	_N = _Nm/_m;
}

vector<tMatrix> KnownChannelOrderAlgorithm::ProcessTrainingSequence(tMatrix observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
// 	int lengthSequenceToProcess = trainingSequence.cols() + _preamble.cols();
	tMatrix toProcessSequence = Util::Append(_preamble,trainingSequence);
	int lengthToProcessSequence = toProcessSequence.cols();
	
	if(observations.cols()<lengthToProcessSequence)
		throw RuntimeException("Insufficient number of observations.");
	
	vector<tMatrix> estimatedMatrices(lengthToProcessSequence);

	// selects all the rows from a symbols matrix
	tRange allSymbolRows(0,_N-1);

	for(int i=_m-1;i<lengthToProcessSequence;i++)
	{
		tRange mColumns(i-_m+1,i);
		estimatedMatrices[i] = _channelEstimator->NextMatrix(observations.col(i),toProcessSequence(allSymbolRows,mColumns),noiseVariances[i]);
	}
	return estimatedMatrices;
}

tMatrix KnownChannelOrderAlgorithm::HsToStackedH(vector<tMatrix> matrices)
{
	
	if((matrices[0].cols() % _m)!=0)
		throw RuntimeException("KnownChannelOrderAlgorithm::HsToStackedH: Incorrect number of columns in the matrices.");
	int L = matrices[0].rows();
	int N = matrices[0].cols()/_m;
	int d = matrices.size()-1;
	
	tMatrix res(matrices[0].rows()*(d+1),N*(_m+d));
    res = 0.0;

	for(int i=0;i<=d;i++)
	{
		tRange rowsRange(i*L,(i+1)*L-1);
		tRange colsRange(i*N,i*N+N*_m-1);
		res(rowsRange,colsRange).inject(matrices[i]);
	}

	return res;
}

