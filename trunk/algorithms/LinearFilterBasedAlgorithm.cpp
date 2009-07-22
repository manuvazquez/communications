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
#include "LinearFilterBasedAlgorithm.h"

// #define DEBUG
#include <KnownSymbolsKalmanEstimator.h>

LinearFilterBasedAlgorithm::LinearFilterBasedAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int backwardsSmoothingLag, int smoothingLag, LinearDetector *linearDetector,  double ARcoefficient, bool substractContributionFromKnownSymbols): KnownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble),_c(backwardsSmoothingLag),_d(smoothingLag),_linearDetector(linearDetector->clone()),_detectedSymbolVectors(N,iLastSymbolVectorToBeDetected),_ARcoefficient(ARcoefficient),_substractContributionFromKnownSymbols(substractContributionFromKnownSymbols)
{
	_estimatedChannelMatrices = new tMatrix[iLastSymbolVectorToBeDetected];
}


LinearFilterBasedAlgorithm::~LinearFilterBasedAlgorithm()
{
	delete _linearDetector;
	delete[] _estimatedChannelMatrices;
}

void LinearFilterBasedAlgorithm::run(tMatrix observations,vector<double> noiseVariances)
{
	process(observations,noiseVariances,tMatrix());
}

void LinearFilterBasedAlgorithm::run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
	process(observations,noiseVariances,trainingSequence);
}

void LinearFilterBasedAlgorithm::process(const tMatrix &observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
	if(observations.rows()!=_nOutputs || trainingSequence.rows()!=_nInputs)
		throw RuntimeException("LinearFilterBasedAlgorithm::process: Observations matrix or training sequence dimensions are wrong.");

	int startDetectionTime = _preamble.cols() + trainingSequence.cols();
	int nObservations = observations.cols();

    if(nObservations<(startDetectionTime+1+_d))
        throw RuntimeException("LinearFilterBasedAlgorithm::process: Not enough observations.");

    vector<tMatrix> trainingSequenceChannelMatrices = _channelEstimator->nextMatricesFromObservationsSequence(observations,noiseVariances,Util::append(_preamble,trainingSequence),_preamble.cols(),startDetectionTime);

    _linearDetector->stateStepsFromObservationsSequence(observations,_d,_preamble.cols(),startDetectionTime);

	for(int j=_preamble.cols();j<startDetectionTime;j++)
	{
		_detectedSymbolVectors.col(j).inject(trainingSequence.col(j-_preamble.cols()));
		_estimatedChannelMatrices[j] = trainingSequenceChannelMatrices[j-_preamble.cols()];
	}

    if(_substractContributionFromKnownSymbols)
    {
        if(_linearDetector->channelMatrixcols() != _nInputs*(_d+1))
            throw RuntimeException("LinearFilterBasedAlgorithm::process: the algorithm is supposed to operate substracting the contribution of the known symbols but this is not compatible with the current linear detector.");
    }

	vector<tMatrix> matricesToStack(_c+_d+1);
	int iSmoothing,iRow;
	tRange rAllObservationsRows(0,_nOutputs-1),rAllSymbolRows(0,_nInputs-1);
	tMatrix stackedNoiseCovariance = LaGenMatDouble::zeros(_nOutputs*(_c+_d+1),_nOutputs*(_c+_d+1));
	double ARcoefficientPower;

	for(int iObservationToBeProcessed=startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
	{
		// already estimated channel matrices are stored in a vector in order to stack them
		for(iSmoothing=-_c;iSmoothing<0;iSmoothing++)
		{
			matricesToStack[iSmoothing+_c] = _estimatedChannelMatrices[iObservationToBeProcessed+iSmoothing];
		}

		ARcoefficientPower = _ARcoefficient;
		for(iSmoothing=0;iSmoothing<=_d;iSmoothing++)
		{
			matricesToStack[_c+iSmoothing] = _estimatedChannelMatrices[iObservationToBeProcessed-1];
			matricesToStack[_c+iSmoothing] *= ARcoefficientPower;
			ARcoefficientPower *= _ARcoefficient;
		}

		tMatrix stackedChannelMatrix = channelMatrices2stackedChannelMatrix(matricesToStack);

		// observation matrix columns that are involved in the smoothing
		tRange rSmoothingRange(iObservationToBeProcessed-_c,iObservationToBeProcessed+_d);

		// the stacked observations vector
		tVector stackedObservations = Util::toVector(observations(rAllObservationsRows,rSmoothingRange),columnwise);

		// stacked noise covariance needs to be constructed
		for(iSmoothing=-_c;iSmoothing<=_d;iSmoothing++)
			for(iRow=0;iRow<_nOutputs;iRow++)
				stackedNoiseCovariance((iSmoothing+_c)*_nOutputs+iRow,(iSmoothing+_c)*_nOutputs+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];

        tVector softEstimations;

        if(_substractContributionFromKnownSymbols)
        {
            softEstimations =  _linearDetector->detect(
                // the last range chooses all the already detected symbol vectors
                substractKnownSymbolsContribution(matricesToStack,_channelOrder,_c,_d,stackedObservations,_detectedSymbolVectors(rAllSymbolRows,tRange(iObservationToBeProcessed-_c-_channelOrder+1,iObservationToBeProcessed-1))),
                // only a part of the channel matrix is needed. The first range chooses all the stacked observation rows
                stackedChannelMatrix(tRange(0,_nOutputs*(_c+_d+1)-1),tRange((_c+_channelOrder-1)*_nInputs,stackedChannelMatrix.cols()-1)),
                stackedNoiseCovariance);
        } else
        {
            softEstimations =  _linearDetector->detect(stackedObservations,stackedChannelMatrix,stackedNoiseCovariance);
        }

		for(iRow=0;iRow<_nInputs;iRow++)
			_detectedSymbolVectors(iRow,iObservationToBeProcessed) = _alphabet.hardDecision(softEstimations(iRow));

		tRange rInvolvedSymbolVectors(iObservationToBeProcessed-_channelOrder+1,iObservationToBeProcessed);
		_estimatedChannelMatrices[iObservationToBeProcessed] = _channelEstimator->nextMatrix(observations.col(iObservationToBeProcessed),_detectedSymbolVectors(rAllSymbolRows,rInvolvedSymbolVectors),noiseVariances[iObservationToBeProcessed]);
	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

tMatrix LinearFilterBasedAlgorithm::getDetectedSymbolVectors()
{
	return _detectedSymbolVectors(tRange(0,_nInputs-1),tRange(_preamble.cols(),_iLastSymbolVectorToBeDetected-1));
}

vector<tMatrix> LinearFilterBasedAlgorithm::getEstimatedChannelMatrices()
{
    vector<tMatrix> channelMatrices;
    channelMatrices.reserve(_iLastSymbolVectorToBeDetected-_preamble.cols());

    for(int i=_preamble.cols();i<_iLastSymbolVectorToBeDetected;i++)
	    channelMatrices.push_back(_estimatedChannelMatrices[i]);

    return channelMatrices;
}
