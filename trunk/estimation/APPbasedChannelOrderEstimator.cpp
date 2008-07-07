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

APPbasedChannelOrderEstimator::APPbasedChannelOrderEstimator(int N,std::vector<int> candidateOrders): ChannelOrderEstimator(N,candidateOrders),_unnormalizedChannelOrderAPPs(candidateOrders.size()),_maxChannelOrder(candidateOrders[Util::Max(candidateOrders)]),_NmaxChannelOrder(_N*_maxChannelOrder),_channelOrder2index(_maxChannelOrder+1,-1),_symbolVector(_NmaxChannelOrder)
{
	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
		_channelOrder2index[_candidateOrders[iChannelOrder]] = iChannelOrder;
}


APPbasedChannelOrderEstimator::~APPbasedChannelOrderEstimator()
{
}

APPbasedChannelOrderEstimator* APPbasedChannelOrderEstimator::Clone()
{
	return new APPbasedChannelOrderEstimator(*this);
}


tMatrix APPbasedChannelOrderEstimator::ComputeProbabilities(const tMatrix& observations,const std::vector<std::vector<tMatrix> > &channelMatrices,const std::vector< double > &noiseVariances,const tMatrix &sequenceToProcess, int iFrom)
{
    int nProbabilitiesToCompute = sequenceToProcess.cols() - iFrom;
    double normalizationCt;
    tMatrix computedChannelOrderAPPs(_candidateOrders.size(),nProbabilitiesToCompute);

    if(observations.cols() < sequenceToProcess.cols())
        throw RuntimeException("APPbasedChannelOrderEstimator::ComputeProbabilities: Insufficient number of observations.");

    if(channelMatrices[0].size() < nProbabilitiesToCompute)
   		throw RuntimeException("APPbasedChannelOrderEstimator::ComputeProbabilities: Insufficient number of channel matrices per channel order.");

    tVector noiselessObservation(observations.rows());

    uint iChannelOrder;
    for(int i=iFrom;i<sequenceToProcess.cols();i++)
    {
        normalizationCt = 0.0;

        _symbolVector = Util::ToVector(sequenceToProcess(_rAllSymbolRows,tRange(i-_maxChannelOrder+1,i)),columnwise);

        for(iChannelOrder=0;iChannelOrder<channelMatrices.size();iChannelOrder++)
        {
            tRange rInvolvedSymbols(_NmaxChannelOrder-_candidateOrders[iChannelOrder]*_N,_NmaxChannelOrder-1);

            // noiselessObservation = LastEstimatedChannelMatrix * stackedSymbolVector
            Blas_Mat_Vec_Mult(channelMatrices[iChannelOrder][i-iFrom],_symbolVector(rInvolvedSymbols),noiselessObservation);

            _unnormalizedChannelOrderAPPs[iChannelOrder] = _channelOrderAPPs[iChannelOrder]* StatUtil::NormalPdf(observations.col(i),noiselessObservation,noiseVariances[i]);

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

void APPbasedChannelOrderEstimator::Update(const tVector &observations,const vector<tMatrix> &channelMatrices,const tVector &symbolVector,double noiseVariance)
{
	if(channelMatrices.size()!=_candidateOrders.size())
		throw RuntimeException("APPbasedChannelOrderEstimator::Update: insufficient number of channel matrices.");

	if(symbolVector.size()!=_N)
		throw RuntimeException("APPbasedChannelOrderEstimator::Update: symbols vector does not have the proper size.");

	// the vector with the symbols involved in the observation is shifted by N
	Util::ShiftUp(_symbolVector,_N);

	// the just detected vector is stored at the end
	for(int i=0;i<_N;i++)
		_symbolVector(_NmaxChannelOrder-_N+i) = symbolVector(i);

	tVector noiselessObservation(observations.size());

    double normalizationCt = 0.0;

	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
	{
		if(channelMatrices[iChannelOrder].cols()!=(_N*_candidateOrders[iChannelOrder]))
			throw RuntimeException("APPbasedChannelOrderEstimator::Update: one (or several) channel matrices has not the proper dimensions.");

		tRange rInvolvedSymbols(_NmaxChannelOrder-_candidateOrders[iChannelOrder]*_N,_NmaxChannelOrder-1);

		Blas_Mat_Vec_Mult(channelMatrices[iChannelOrder],_symbolVector(rInvolvedSymbols),noiselessObservation);

		_unnormalizedChannelOrderAPPs[iChannelOrder] = _channelOrderAPPs[iChannelOrder]* StatUtil::NormalPdf(observations,noiselessObservation,noiseVariance);

		normalizationCt += _unnormalizedChannelOrderAPPs[iChannelOrder];
	}

	if(normalizationCt!=0.0)
		for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
		{
			_channelOrderAPPs[iChannelOrder] = _unnormalizedChannelOrderAPPs[iChannelOrder] / normalizationCt;
		}
}
