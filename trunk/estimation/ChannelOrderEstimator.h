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
#ifndef CHANNELORDERESTIMATOR_H
#define CHANNELORDERESTIMATOR_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <vector>

class ChannelOrderEstimator{
protected:
	int _N;
	tMatrix _preamble;
	std::vector<int> _candidateOrders;
    std::vector<double> _channelOrderAPPs;
public:
    ChannelOrderEstimator(int N, const tMatrix &preamble, std::vector<int> candidateOrders);

    ChannelOrderEstimator(const tMatrix &preamble, std::vector<int> candidateOrders, std::vector<double> channelOrderAPPs);

    virtual ~ChannelOrderEstimator() {}

    double GetChannelOrderAPP(int n) {return _channelOrderAPPs[n];}

    tVector GetChannelOrderAPPsVector();

    virtual ChannelOrderEstimator *Clone() = 0;

	virtual void Update(const tVector &observations,const std::vector<tMatrix> &channelMatrix,const tVector &symbolsVector,double noiseVariance)=0;

	virtual std::vector<double> ComputeProbabilities(const tMatrix& observations,const std::vector<std::vector<tMatrix> > channelMatrices, std::vector< double > noiseVariances, tMatrix symbolVectors)=0;

};

#endif
