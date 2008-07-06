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
#ifndef APPBASEDCHANNELORDERESTIMATOR_H
#define APPBASEDCHANNELORDERESTIMATOR_H

#include <ChannelOrderEstimator.h>

/**
    @author Manu <manu@rustneversleeps>
*/

#include <ChannelMatrixEstimator.h>
#include <Util.h>
#include <StatUtil.h>
#include <lapackpp/blas2pp.h>

class APPbasedChannelOrderEstimator : public ChannelOrderEstimator
{
protected:
    tRange _rAllSymbolRows;
    vector<double> _unnormalizedChannelOrderAPPs;
    int _maxChannelOrder,_NmaxChannelOrder;
    vector<int> _channelOrder2index;
    tVector _symbolVector;
public:
    APPbasedChannelOrderEstimator(int N, const tMatrix& preamble, std::vector<int> candidateOrders);

    ~APPbasedChannelOrderEstimator();

    virtual APPbasedChannelOrderEstimator *Clone();

    void Update(const tVector &observations,const vector<tMatrix> &channelMatrix,const tVector &symbolVector,double noiseVariance);

    virtual tMatrix ComputeProbabilities(const tMatrix& observations,const std::vector<std::vector<tMatrix> > channelMatrices, std::vector< double > noiseVariances, tMatrix sequenceToProcess, int iFrom);

};

#endif
