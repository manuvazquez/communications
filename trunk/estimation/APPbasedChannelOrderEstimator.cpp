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

APPbasedChannelOrderEstimator::APPbasedChannelOrderEstimator(const tMatrix& preamble, std::vector<int> candidateOrders, vector<tMatrix> initialChannelMatrixEstimations,double ARcoefficient): ChannelOrderEstimator(preamble,candidateOrders),_lastEstimatedChannelMatrices(initialChannelMatrixEstimations),_rAllSymbolRows(0,_preamble.rows()-1),_unnormalizedChannelOrderAPPs(initialChannelMatrixEstimations.size()),_ARcoefficient(ARcoefficient)
{
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

        for(iChannelOrder=0;iChannelOrder<channelMatrices.size();iChannelOrder++)
        {
            tRange rInvolvedSymbolVectors(i-_candidateOrders[iChannelOrder]+1,i);
            tVector stackedSymbolVector = Util::ToVector(sequenceToProcess(_rAllSymbolRows,rInvolvedSymbolVectors),columnwise);

            tMatrix predictedChannelMatrix = _lastEstimatedChannelMatrices[iChannelOrder];
            predictedChannelMatrix *= _ARcoefficient;

			#ifdef DEBUG
				cout << "Antes de llamar a Blas" << endl;
			#endif

            // predictedNoiselessObservation = LastEstimatedChannelMatrix * stackedSymbolVector
            Blas_Mat_Vec_Mult(predictedChannelMatrix,stackedSymbolVector,predictedNoiselessObservation);

            _unnormalizedChannelOrderAPPs[iChannelOrder] = _channelOrderAPPs[iChannelOrder]* StatUtil::NormalPdf(observations.col(i),predictedNoiselessObservation,noiseVariances[i]);

            normalizationCt += _unnormalizedChannelOrderAPPs[iChannelOrder];

			// the next estimated channel matrix is stored into
            _lastEstimatedChannelMatrices[iChannelOrder] = channelMatrices[iChannelOrder][i-_preamble.cols()];

			#ifdef DEBUG
				cout << "Final del bucle iChannelOrder" << endl;
			#endif
        }

        if(normalizationCt!=0)
            for(uint iChannelOrder=0;iChannelOrder<channelMatrices.size();iChannelOrder++)
            {
                _channelOrderAPPs[iChannelOrder] = _unnormalizedChannelOrderAPPs[iChannelOrder] / normalizationCt;
            }

		#ifdef DEBUG
			cout << "Final del bucle i" << endl;
		#endif
    }

	#ifdef DEBUG
    	cout << "Final del metodo" << endl;
   	#endif

   	return _channelOrderAPPs;
}

