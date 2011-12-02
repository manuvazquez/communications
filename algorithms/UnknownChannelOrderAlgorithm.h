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
#ifndef UNKNOWNCHANNELORDERALGORITHM_H
#define UNKNOWNCHANNELORDERALGORITHM_H

#include <UnknownChannelAlgorithm.h>

/**
    @author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <ChannelMatrixEstimator.h>

class UnknownChannelOrderAlgorithm : public UnknownChannelAlgorithm
{
protected:
    vector<ChannelMatrixEstimator *> _channelEstimators;
    vector<uint> _candidateOrders;
    uint _maxOrder,_iFirstObservation,_nInputsXmaxChannelOrder;
    uint *_channelOrder2index;
    MatrixXd _preamble;   

	// a vector of matrices with the channel order a posteriori probabilities. Each matrix corresponds to a different output (or receiving antenna)
	// if all the outputs have the same channel order (and channel order APPs, consequently), only one is stored, i.e., it's a vector with a single
	// element
    std::vector<MatrixXd> _channelOrderAPPs;
// 	MatrixXd _channelOrderAPPs;
public:
    UnknownChannelOrderAlgorithm(string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected,vector<ChannelMatrixEstimator *> channelEstimators,MatrixXd preamble,uint iFirstObservation);

    ~UnknownChannelOrderAlgorithm();

	// if no output is specified, the first one is assumed
	/*!
	  It returns a matrix that contains the a posteriori channel order probabilities at all time instants for the FIRST output
	  \return matrix in which each row represents a candidate channel order, and each column a time instant
	*/
    MatrixXd getComputedChannelOrderAPPs() { return _channelOrderAPPs[0].block(0,_preamble.cols(),_candidateOrders.size(),_iLastSymbolVectorToBeDetected-_preamble.cols());}

	/*!
	  It returns a matrix that contains the a posteriori channel order probabilities at all time instants for the a given output
	  \param iOutput the output for which the channel order APPs are wanted
	  \return matrix in which each row represents a candidate channel order, and each column a time instant
	*/
    MatrixXd getComputedChannelOrderAPPs(uint iOutput) { return _channelOrderAPPs[iOutput].block(0,_preamble.cols(),_candidateOrders.size(),_iLastSymbolVectorToBeDetected-_preamble.cols());}
    
    virtual bool estimatesOneSingleChannelOrder() const { return true;}
};

#endif
