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
#ifndef KNOWNCHANNELORDERALGORITHM_H
#define KNOWNCHANNELORDERALGORITHM_H

#include <UnknownChannelAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <types.h>
#include <Util.h>

class KnownChannelOrderAlgorithm : public UnknownChannelAlgorithm
{
protected:
	ChannelMatrixEstimator *_channelEstimator;
	int _channelOrder,_nInputsXchannelOrder;
    MatrixXd _preamble;   

public:
    KnownChannelOrderAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m, ChannelMatrixEstimator *channelEstimator,MatrixXd preamble);
    KnownChannelOrderAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m,MatrixXd preamble);
	~KnownChannelOrderAlgorithm();

	using Algorithm::channelMatrices2stackedChannelMatrix;
    MatrixXd channelMatrices2stackedChannelMatrix(vector<MatrixXd> matrices) { return channelMatrices2stackedChannelMatrix(matrices,_channelOrder);}   
};

#endif