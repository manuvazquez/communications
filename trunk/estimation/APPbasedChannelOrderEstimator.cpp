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

APPbasedChannelOrderEstimator::APPbasedChannelOrderEstimator(int N, const tMatrix& preamble, std::vector<int> candidateOrders, vector<tMatrix> initialChannelMatrixEstimations,double ARcoefficient): ChannelOrderEstimator(N,preamble,candidateOrders),_lastEstimatedChannelMatrices(initialChannelMatrixEstimations),_rAllSymbolRows(0,_preamble.rows()-1),_unnormalizedChannelOrderAPPs(initialChannelMatrixEstimations.size()),_ARcoefficient(ARcoefficient),_maxChannelOrder(candidateOrders[Util::Max(candidateOrders)]),_NmaxChannelOrder(_N*_maxChannelOrder),_channelOrder2index(_maxChannelOrder+1,-1),_symbolsVector(_NmaxChannelOrder)
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


vector< double > APPbasedChannelOrderEstimator::ComputeProbabilities(const tMatrix& observations,const vector<vector<tMatrix> > channelMatrices, vector< double > noiseVariances, tMatrix symbolVectors)
{
    tMatrix sequenceToProcess = Util::Append(_preamble,symbolVectors);
    int lengthSequenceToProcess = sequenceToProcess.cols();
    double normalizationCt;

    if(observations.cols() < lengthSequenceToProcess)
        throw RuntimeException("APPbasedChannelOrderEstimator::ComputeProbabilities: Insufficient number of observations.");

	if(channelMatrices.size() != _lastEstimatedChannelMatrices.size())
   		throw RuntimeException("APPbasedChannelOrderEstimator::ComputeProbabilities: \"channelMatrices\" size is not coherent with received initial channel matrix estimations.");

    if(channelMatrices[0].size() < symbolVectors.cols())
   		throw RuntimeException("APPbasedChannelOrderEstimator::ComputeProbabilities: Insufficient number of channel matrices per channel order.");

    tVector predictedNoiselessObservation(observations.rows());

    uint iChannelOrder;
    for(int i=_preamble.cols();i<lengthSequenceToProcess;i++)
    {
        normalizationCt = 0.0;

        _symbolsVector = Util::ToVector(sequenceToProcess(_rAllSymbolRows,tRange(i-_maxChannelOrder+1,i)),columnwise);

        for(iChannelOrder=0;iChannelOrder<channelMatrices.size();iChannelOrder++)
        {
            tRange rInvolvedSymbols(_NmaxChannelOrder-_candidateOrders[iChannelOrder]*_N,_NmaxChannelOrder-1);

            tMatrix predictedChannelMatrix = _lastEstimatedChannelMatrices[iChannelOrder];
            predictedChannelMatrix *= _ARcoefficient;

            // predictedNoiselessObservation = LastEstimatedChannelMatrix * stackedSymbolVector
            Blas_Mat_Vec_Mult(predictedChannelMatrix,_symbolsVector(rInvolvedSymbols),predictedNoiselessObservation);

            _unnormalizedChannelOrderAPPs[iChannelOrder] = _channelOrderAPPs[iChannelOrder]* StatUtil::NormalPdf(observations.col(i),predictedNoiselessObservation,noiseVariances[i]);

            normalizationCt += _unnormalizedChannelOrderAPPs[iChannelOrder];

			// the next estimated channel matrix is stored into
            _lastEstimatedChannelMatrices[iChannelOrder] = channelMatrices[iChannelOrder][i-_preamble.cols()];
        }

        if(normalizationCt!=0.0)
            for(uint iChannelOrder=0;iChannelOrder<channelMatrices.size();iChannelOrder++)
            {
                _channelOrderAPPs[iChannelOrder] = _unnormalizedChannelOrderAPPs[iChannelOrder] / normalizationCt;
            }
    }

   	return _channelOrderAPPs;
}

void APPbasedChannelOrderEstimator::Update(const tVector &observations,const vector<tMatrix> &channelMatrices,const tVector &symbolsVector,double noiseVariance)
{
	if(channelMatrices.size()!=_candidateOrders.size())
		throw RuntimeException("APPbasedChannelOrderEstimator::Update: insufficient number of channel matrices.");
// 	if(channelMatrix.cols()!=symbolsVector.size())
// 		throw RuntimeException("APPbasedChannelOrderEstimator::Update: channel matrix and symbols vector dimensions are not coherent.");

// 	int iChannelOrder = _channelOrder2index[channelMatrix.cols()/_preamble.rows()];
// 	if(iChannelOrder==-1)
// 		throw RuntimeException("APPbasedChannelOrderEstimator::Update: the computed order of the channel matrix is none of the candidates.");

	if(symbolsVector.size()!=_N)
		throw RuntimeException("APPbasedChannelOrderEstimator::Update: symbols vector does not have the proper size.");

	// the vector with the symbols involved in the observation is shifted by N
	Util::ShiftUp(_symbolsVector,_N);

	// the just detected vector is stored at the end
	for(int i=0;i<_N;i++)
		_symbolsVector(_NmaxChannelOrder-_N+i) = symbolsVector(i);

	tVector predictedNoiselessObservation(observations.size());

    double normalizationCt = 0.0;

	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
	{
		if(channelMatrices[iChannelOrder].cols()!=(_N*_candidateOrders[iChannelOrder]))
			throw RuntimeException("APPbasedChannelOrderEstimator::Update: one (or several) channel matrices has not the proper dimensions.");

		tRange rInvolvedSymbols(_NmaxChannelOrder-_candidateOrders[iChannelOrder]*_N,_NmaxChannelOrder-1);

		Blas_Mat_Vec_Mult(channelMatrices[iChannelOrder],_symbolsVector(rInvolvedSymbols),predictedNoiselessObservation);

		_unnormalizedChannelOrderAPPs[iChannelOrder] = _channelOrderAPPs[iChannelOrder]* StatUtil::NormalPdf(observations,predictedNoiselessObservation,noiseVariance);

		normalizationCt += _unnormalizedChannelOrderAPPs[iChannelOrder];
	}

	if(normalizationCt!=0.0)
		for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
		{
			_channelOrderAPPs[iChannelOrder] = _unnormalizedChannelOrderAPPs[iChannelOrder] / normalizationCt;
		}
}
