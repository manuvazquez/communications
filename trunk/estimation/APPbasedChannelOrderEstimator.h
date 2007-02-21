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
    vector<tMatrix> _lastEstimatedChannelMatrices;
    tRange _rAllSymbolRows;
    vector<double> _unnormalizedChannelOrderAPPs;
    double _ARcoefficient;
public:
    APPbasedChannelOrderEstimator(const tMatrix& preamble, std::vector<int> candidateOrders, vector<tMatrix> initialChannelMatrixEstimations,double ARcoefficient);

    ~APPbasedChannelOrderEstimator();

    virtual APPbasedChannelOrderEstimator *Clone();

    virtual std::vector< double > ComputeProbabilities(const tMatrix& observations, const vector<vector<tMatrix> > channelMatrices, std::vector< double > noiseVariances, tMatrix symbolVectors);

};

#endif
