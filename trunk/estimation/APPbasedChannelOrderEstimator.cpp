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

APPbasedChannelOrderEstimator::APPbasedChannelOrderEstimator(const tMatrix& preamble,vector<ChannelMatrixEstimator *> channelEstimators,double ARcoefficient): ChannelOrderEstimator(preamble),_channelEstimators(channelEstimators.size()),_rAllSymbolRows(0,_preamble.rows()-1),_channelOrderAPPs(channelEstimators.size(),1.0/(double)channelEstimators.size()),_unnormalizedChannelOrderAPPs(new double[channelEstimators.size()]),_ARcoefficient(ARcoefficient)
{
    for(int i=0;i<channelEstimators.size();i++)
        _channelEstimators[i] = channelEstimators[i]->Clone();
}


APPbasedChannelOrderEstimator::~APPbasedChannelOrderEstimator()
{
    for(int i=0;i<_channelEstimators.size();i++)
        delete _channelEstimators[i];

    delete[] _unnormalizedChannelOrderAPPs;
}


vector< double > APPbasedChannelOrderEstimator::ComputeProbabilities(const tMatrix& observations, vector< double > noiseVariances, tMatrix symbolVectors)
{
    tMatrix sequenceToProcess = Util::Append(_preamble,symbolVectors);
    int lengthSequenceToProcess = sequenceToProcess.cols();
    double normalizationCt;

    if(observations.cols() < lengthSequenceToProcess)
        throw RuntimeException("APPbasedChannelOrderEstimator::ComputeProbabilities: Insufficient number of observations.");

    tVector predictedNoiselessObservation(observations.rows());

    int iChannelOrder;
    for(int i=_preamble.cols();i<lengthSequenceToProcess;i++)
    {
        normalizationCt = 0.0;

        for(iChannelOrder=0;iChannelOrder<_channelEstimators.size();iChannelOrder++)
        {
            tRange rInvolvedSymbolVectors(i-_channelEstimators[iChannelOrder]->Memory()+1,i);
            tVector stackedSymbolVector = Util::ToVector(sequenceToProcess(_rAllSymbolRows,rInvolvedSymbolVectors),columnwise);

            tMatrix predictedChannelMatrix = _channelEstimators[iChannelOrder]->LastEstimatedChannelMatrix();
            predictedChannelMatrix *= _ARcoefficient;

            // predictedNoiselessObservation = LastEstimatedChannelMatrix * stackedSymbolVector
            Blas_Mat_Vec_Mult(predictedChannelMatrix,stackedSymbolVector,predictedNoiselessObservation);

            _unnormalizedChannelOrderAPPs[iChannelOrder] = _channelOrderAPPs[iChannelOrder]* StatUtil::NormalPdf(observations.col(i),predictedNoiselessObservation,noiseVariances[i]);

            normalizationCt += _unnormalizedChannelOrderAPPs[iChannelOrder];

            _channelEstimators[iChannelOrder]->NextMatrix(observations.col(i),sequenceToProcess(_rAllSymbolRows,rInvolvedSymbolVectors),noiseVariances[i]);
        }

        if(normalizationCt!=0)
            for(int iChannelOrder=0;iChannelOrder<_channelEstimators.size();iChannelOrder++)
            {
                _channelOrderAPPs[iChannelOrder] = _unnormalizedChannelOrderAPPs[iChannelOrder] / normalizationCt;
            }
    }
}

