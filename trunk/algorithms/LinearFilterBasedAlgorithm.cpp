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

LinearFilterBasedAlgorithm::LinearFilterBasedAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int backwardsSmoothingLag, int smoothingLag, LinearDetector *linearDetector,  double ARcoefficient, bool substractContributionFromKnownSymbols): KnownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble),_c(backwardsSmoothingLag),_d(smoothingLag),_linearDetector(linearDetector->clone()),_detectedSymbolVectors(N,iLastSymbolVectorToBeDetected),_estimatedChannelMatrices(iLastSymbolVectorToBeDetected),_ARcoefficient(ARcoefficient),_substractContributionFromKnownSymbols(substractContributionFromKnownSymbols)
{
}


LinearFilterBasedAlgorithm::~LinearFilterBasedAlgorithm()
{
	delete _linearDetector;
}

void LinearFilterBasedAlgorithm::run(MatrixXd observations,vector<double> noiseVariances)
{
    process(observations,noiseVariances,MatrixXd());
}

void LinearFilterBasedAlgorithm::run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence)
{
    process(observations,noiseVariances,trainingSequence);
}

void LinearFilterBasedAlgorithm::process(const MatrixXd &observations,vector<double> noiseVariances, MatrixXd trainingSequence)
{
    if(observations.rows()!=_nOutputs || trainingSequence.rows()!=_nInputs)
        throw RuntimeException("LinearFilterBasedAlgorithm::process: Observations matrix or training sequence dimensions are wrong.");

    int startDetectionTime = _preamble.cols() + trainingSequence.cols();
    int nObservations = observations.cols();

    if(nObservations<(startDetectionTime+1+_d))
        throw RuntimeException("LinearFilterBasedAlgorithm::process: Not enough observations.");

    MatrixXd preambleTrainingSequence(_preamble.rows(),_preamble.cols()+trainingSequence.cols());
    preambleTrainingSequence << Util::lapack2eigen(_preamble),trainingSequence;

    vector<MatrixXd> trainingSequenceChannelMatrices = _channelEstimator->nextMatricesFromObservationsSequence(observations,noiseVariances,preambleTrainingSequence,_preamble.cols(),startDetectionTime);

    _linearDetector->stateStepsFromObservationsSequence(observations,_d,_preamble.cols(),startDetectionTime);

    for(int j=_preamble.cols();j<startDetectionTime;j++)
    {
        _detectedSymbolVectors.col(j) = trainingSequence.col(j-_preamble.cols());
        _estimatedChannelMatrices[j] = trainingSequenceChannelMatrices[j-_preamble.cols()];
    }

    if(_substractContributionFromKnownSymbols)
    {
        if(_linearDetector->channelMatrixcols() != _nInputs*(_d+1))
            throw RuntimeException("LinearFilterBasedAlgorithm::process: the algorithm is supposed to operate substracting the contribution of the known symbols but this is not compatible with the current linear detector.");
    }

    vector<MatrixXd> matricesToStack(_c+_d+1);
    int iSmoothing,iRow;
    MatrixXd stackedNoiseCovariance = MatrixXd::Zero(_nOutputs*(_c+_d+1),_nOutputs*(_c+_d+1));
    double ARcoefficientPower;

    for(int iObservationToBeProcessed=startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
        // already estimated channel matrices are stored in a vector in order to stack them
        for(iSmoothing=-_c;iSmoothing<0;iSmoothing++)
            matricesToStack[iSmoothing+_c] = _estimatedChannelMatrices[iObservationToBeProcessed+iSmoothing];

        ARcoefficientPower = _ARcoefficient;
        for(iSmoothing=0;iSmoothing<=_d;iSmoothing++)
        {
            matricesToStack[_c+iSmoothing] = _estimatedChannelMatrices[iObservationToBeProcessed-1];
            matricesToStack[_c+iSmoothing] *= ARcoefficientPower;
            ARcoefficientPower *= _ARcoefficient;
        }

        MatrixXd stackedChannelMatrix = channelMatrices2stackedChannelMatrix(matricesToStack);

        VectorXd stackedObservations = Util::toVector(observations.block(0,iObservationToBeProcessed-_c,_nOutputs,_c+_d+1),columnwise);

        // stacked noise covariance needs to be constructed
        for(iSmoothing=-_c;iSmoothing<=_d;iSmoothing++)
            for(iRow=0;iRow<_nOutputs;iRow++)
                stackedNoiseCovariance((iSmoothing+_c)*_nOutputs+iRow,(iSmoothing+_c)*_nOutputs+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];

        VectorXd softEstimations;

        if(_substractContributionFromKnownSymbols)
        {
            softEstimations =  _linearDetector->detect(
                // the last range chooses all the already detected symbol vectors
                substractKnownSymbolsContribution(matricesToStack,_channelOrder,_c,_d,stackedObservations,_detectedSymbolVectors.block(0,iObservationToBeProcessed-_c-_channelOrder+1,_nInputs,_c+_channelOrder-1)),
                // only a part of the channel matrix is needed. The first range chooses all the stacked observation rows
                stackedChannelMatrix.block(0,(_c+_channelOrder-1)*_nInputs,_nOutputs*(_c+_d+1),stackedChannelMatrix.cols()-(_c+_channelOrder-1)*_nInputs),
                stackedNoiseCovariance);                
        } else
            softEstimations =  _linearDetector->detect(stackedObservations,stackedChannelMatrix,stackedNoiseCovariance);

        for(iRow=0;iRow<_nInputs;iRow++)
            _detectedSymbolVectors(iRow,iObservationToBeProcessed) = _alphabet.hardDecision(softEstimations(iRow));

        _estimatedChannelMatrices[iObservationToBeProcessed] = _channelEstimator->nextMatrix(observations.col(iObservationToBeProcessed),_detectedSymbolVectors.block(0,iObservationToBeProcessed-_channelOrder+1,_nInputs,_channelOrder),noiseVariances[iObservationToBeProcessed]);
    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

MatrixXd LinearFilterBasedAlgorithm::getDetectedSymbolVectors_eigen()
{
    return _detectedSymbolVectors.block(0,_preamble.cols(),_nInputs,_iLastSymbolVectorToBeDetected-_preamble.cols());
}

vector<MatrixXd> LinearFilterBasedAlgorithm::getEstimatedChannelMatrices_eigen()
{
    vector<MatrixXd> channelMatrices;
    channelMatrices.reserve(_iLastSymbolVectorToBeDetected-_preamble.cols());

    for(int i=_preamble.cols();i<_iLastSymbolVectorToBeDetected;i++)
        channelMatrices.push_back(_estimatedChannelMatrices[i]);

    return channelMatrices;
}
