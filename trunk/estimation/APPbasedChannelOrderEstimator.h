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
    vector<ChannelMatrixEstimator *> _channelEstimators;
    vector<double> _channelOrderAPPs;
    tRange _rAllSymbolRows;
    double *_unnormalizedChannelOrderAPPs;
    double _ARcoefficient;
public:
    APPbasedChannelOrderEstimator(const tMatrix& preamble,vector<ChannelMatrixEstimator *> channelEstimators,double ARcoefficient);

    ~APPbasedChannelOrderEstimator();

    virtual std::vector< double > ComputeProbabilities(const tMatrix& observations, std::vector< double > noiseVariances, tMatrix symbolVectors);

};

#endif
