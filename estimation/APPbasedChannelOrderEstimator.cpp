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

// #define DEBUG

APPbasedChannelOrderEstimator::APPbasedChannelOrderEstimator(uint N,std::vector<uint> candidateOrders): ChannelOrderEstimator(N,candidateOrders),_unnormalizedChannelOrderAPPs(candidateOrders.size()),_maxChannelOrder(candidateOrders[Util::max(candidateOrders)]),_nInputs_maxChannelOrder(_nInputs*_maxChannelOrder),_channelOrder2index(_maxChannelOrder+1,-1)/*,_symbolVector(_nInputs_maxChannelOrder)*/
{
    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _channelOrder2index[_candidateOrders[iChannelOrder]] = iChannelOrder;
}

APPbasedChannelOrderEstimator* APPbasedChannelOrderEstimator::clone()
{
    return new APPbasedChannelOrderEstimator(*this);
}

MatrixXd APPbasedChannelOrderEstimator::computeProbabilities(const MatrixXd& observations,const std::vector<std::vector<MatrixXd> > &channelMatrices,const std::vector< double > &noiseVariances,const MatrixXd &sequenceToProcess, int iFrom)
{
    uint nProbabilitiesToCompute = sequenceToProcess.cols() - iFrom;
    double normalizationCt;
    MatrixXd computedChannelOrderAPPs = MatrixXd::Zero(_candidateOrders.size(),nProbabilitiesToCompute);

    if(observations.cols() < sequenceToProcess.cols())
        throw RuntimeException("APPbasedChannelOrderEstimator::computeProbabilities: insufficient number of observations.");

    if(channelMatrices[0].size() < nProbabilitiesToCompute)
        throw RuntimeException("APPbasedChannelOrderEstimator::computeProbabilities: insufficient number of channel matrices per channel order.");

    uint iChannelOrder;
    for(int i=iFrom;i<sequenceToProcess.cols();i++)
    {
        normalizationCt = 0.0;

        _symbolVector = Util::toVector(sequenceToProcess.block(0,i-_maxChannelOrder+1,_nInputs,_maxChannelOrder),columnwise);

        for(iChannelOrder=0;iChannelOrder<channelMatrices.size();iChannelOrder++)
        {
            _unnormalizedChannelOrderAPPs[iChannelOrder] = _channelOrderAPPs[iChannelOrder]* StatUtil::normalPdf(observations.col(i),channelMatrices[iChannelOrder][i-iFrom]*_symbolVector.tail(_candidateOrders[iChannelOrder]*_nInputs),noiseVariances[i]);

            normalizationCt += _unnormalizedChannelOrderAPPs[iChannelOrder];
        }

        if(normalizationCt!=0.0)
            for(uint iChannelOrder=0;iChannelOrder<channelMatrices.size();iChannelOrder++)
            {
                _channelOrderAPPs[iChannelOrder] = _unnormalizedChannelOrderAPPs[iChannelOrder] / normalizationCt;
                computedChannelOrderAPPs(iChannelOrder,i-iFrom) = _channelOrderAPPs[iChannelOrder];
            }
    }

    return computedChannelOrderAPPs;
}

void APPbasedChannelOrderEstimator::update(const VectorXd &observations,const std::vector<MatrixXd> &channelMatrices,const VectorXd &symbolVector,double noiseVariance)
{
    if(channelMatrices.size()!=_candidateOrders.size())
        throw RuntimeException("APPbasedChannelOrderEstimator::Update: insufficient number of channel matrices.");

    if(symbolVector.size()!=_nInputs)
        throw RuntimeException("APPbasedChannelOrderEstimator::Update: symbols vector does not have the proper size.");

    // the vector with the symbols involved in the observation is shifted by N
    Util::shiftUp(_symbolVector,_nInputs);

    // the just detected vector is stored at the end
    for(uint i=0;i<_nInputs;i++)
        _symbolVector(_nInputs_maxChannelOrder-_nInputs+i) = symbolVector(i);

    double normalizationCt = 0.0;

    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
    {
        if(channelMatrices[iChannelOrder].cols()!=(_nInputs*_candidateOrders[iChannelOrder]))
            throw RuntimeException("APPbasedChannelOrderEstimator::update: one (or several) channel matrices has not the proper dimensions.");
        
        _unnormalizedChannelOrderAPPs[iChannelOrder] = _channelOrderAPPs[iChannelOrder]* StatUtil::normalPdf(observations,channelMatrices[iChannelOrder]*_symbolVector.tail(_candidateOrders[iChannelOrder]*_nInputs),noiseVariance);

        normalizationCt += _unnormalizedChannelOrderAPPs[iChannelOrder];
    }

    if(normalizationCt!=0.0)
        for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        {
            _channelOrderAPPs[iChannelOrder] = _unnormalizedChannelOrderAPPs[iChannelOrder] / normalizationCt;
        }
}
