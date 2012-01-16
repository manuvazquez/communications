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
#include <Util.h>

class ChannelOrderEstimator{
protected:
    uint _nInputs;
    std::vector<uint> _candidateOrders;
    std::vector<double> _channelOrderAPPs;
public:
    ChannelOrderEstimator(uint N, std::vector<uint> candidateOrders);

    ChannelOrderEstimator(std::vector<uint> candidateOrders, std::vector<double> channelOrderAPPs);

    virtual ~ChannelOrderEstimator() {}

    double getChannelOrderAPP(uint n) {return _channelOrderAPPs[n];}

    VectorXd getChannelOrderAPPsVector();

    virtual ChannelOrderEstimator *clone() = 0;

    virtual void update(const VectorXd &observations,const std::vector<MatrixXd> &channelMatrix,const VectorXd &symbolVector,double noiseVariance) = 0;

    virtual MatrixXd computeProbabilities(const MatrixXd& observations,const std::vector<std::vector<MatrixXd> > &channelMatrices,const std::vector< double > &noiseVariances,const MatrixXd &sequenceToProcess, uint iFrom) = 0;

};

#endif
