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
// #include <realData.h>

#include <KnownSymbolsKalmanEstimator.h>

LinearFilterBasedAlgorithm::LinearFilterBasedAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, LinearDetector *linearDetector,  std::vector<double> ARcoefficients, bool substractContributionFromKnownSymbols): 
KnownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble),_d(smoothingLag),_linearDetector(linearDetector->clone()),_detectedSymbolVectors(N,iLastSymbolVectorToBeDetected),_estimatedChannelMatrices(iLastSymbolVectorToBeDetected),_ARcoefficients(ARcoefficients),_substractContributionFromKnownSymbols(substractContributionFromKnownSymbols)
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
	// if there is an actual training sequence...
	if(trainingSequence.rows()!=0)
		// we make sure it has the proper number of rows
		assert(trainingSequence.rows()==_nInputs);

	uint startDetectionTime = _preamble.cols() + trainingSequence.cols();

    MatrixXd preambleTrainingSequence(_preamble.rows(),_preamble.cols()+trainingSequence.cols());
	preambleTrainingSequence << _preamble,trainingSequence;

	// several channel estimation steps are taken in a row
    vector<MatrixXd> trainingSequenceChannelMatrices = _channelEstimator->nextMatricesFromObservationsSequence(observations,noiseVariances,preambleTrainingSequence,_preamble.cols(),startDetectionTime);

	// the linear detector is updated accordingly using the same observations
    _linearDetector->stateStepsFromObservationsSequence(observations,_d,_preamble.cols(),startDetectionTime);

	// if there is a training sequence
    for(uint j=_preamble.cols();j<startDetectionTime;j++)
    {
		// the known symbols are stored as detected symbols
        _detectedSymbolVectors.col(j) = trainingSequence.col(j-_preamble.cols());
		
		// the channel estimates obtained during the training sequence (if present) are stored
        _estimatedChannelMatrices[j] = trainingSequenceChannelMatrices[j-_preamble.cols()];
    }
    
    process(observations,noiseVariances,trainingSequence);
}

void LinearFilterBasedAlgorithm::process(const MatrixXd &observations,vector<double> noiseVariances, MatrixXd trainingSequence)
{
	assert(observations.rows()==_nOutputs); 
	
    uint nObservations = observations.cols();

    uint startDetectionTime = _preamble.cols() + trainingSequence.cols();

	assert(nObservations>=startDetectionTime+1+_d);

    if(_substractContributionFromKnownSymbols)
		// the algorithm is supposed to operate substracting the contribution of the known symbols but this is not compatible with the current linear detector
		assert(_linearDetector->channelMatrixcols()==_nInputs*(_d+1));

    uint iSmoothing,iRow;
    MatrixXd stackedNoiseCovariance = MatrixXd::Zero(_nOutputs*(_d+1),_nOutputs*(_d+1));
	
	std::vector<MatrixXd>::reverse_iterator iterMatrices;	
	std::vector<MatrixXd> ARmatricesBuffer(_ARcoefficients.size(),_channelEstimator->lastEstimatedChannelMatrix());

	// the channel matrices estimated during the training sequence are copied into ARmatricesBuffer (last matrix in the vector is the more recent)
	iterMatrices = ARmatricesBuffer.rbegin();
	for(uint iTrainingSeq=0;iTrainingSeq<trainingSequence.cols() && iterMatrices!=ARmatricesBuffer.rend();iTrainingSeq++,iterMatrices++)
		*iterMatrices = _estimatedChannelMatrices[startDetectionTime-1-iTrainingSeq];
	
    for(uint iObservationToBeProcessed=startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
		vector<MatrixXd> matricesToStack = getChannelMatricesToStackForSmoothing(ARmatricesBuffer);
        
        MatrixXd stackedChannelMatrix = channelMatrices2stackedChannelMatrix(matricesToStack);

        VectorXd stackedObservations = Util::toVector(observations.block(0,iObservationToBeProcessed,_nOutputs,_d+1),columnwise);

        // stacked noise covariance needs to be constructed
        for(iSmoothing=0;iSmoothing<=_d;iSmoothing++)
            for(iRow=0;iRow<_nOutputs;iRow++)
                stackedNoiseCovariance((iSmoothing)*_nOutputs+iRow,(iSmoothing)*_nOutputs+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];

        VectorXd softEstimations;

        if(_substractContributionFromKnownSymbols)
        {
            softEstimations =  _linearDetector->detect(
                // the last range chooses all the already detected symbol vectors
                substractKnownSymbolsContribution(matricesToStack,_channelOrder,_d,stackedObservations,_detectedSymbolVectors.block(0,iObservationToBeProcessed-_channelOrder+1,_nInputs,_channelOrder-1)),
                // only a part of the channel matrix is needed. The first range chooses all the stacked observation rows
                stackedChannelMatrix.block(0,(_channelOrder-1)*_nInputs,_nOutputs*(_d+1),stackedChannelMatrix.cols()-(_channelOrder-1)*_nInputs),
                stackedNoiseCovariance);                
        } else
            softEstimations =  _linearDetector->detect(stackedObservations,stackedChannelMatrix,stackedNoiseCovariance);

        for(iRow=0;iRow<_nInputs;iRow++)
            _detectedSymbolVectors(iRow,iObservationToBeProcessed) = _alphabet.hardDecision(softEstimations(iRow));

		_estimatedChannelMatrices[iObservationToBeProcessed] = _channelEstimator->nextMatrix(observations.col(iObservationToBeProcessed),
																							 obtainChannelMatrixEstimatorFeed(softEstimations,_detectedSymbolVectors.block(0,iObservationToBeProcessed-_channelOrder+1,_nInputs,_channelOrder)),
																							 noiseVariances[iObservationToBeProcessed]);
// 																							 noiseVariances[iObservationToBeProcessed]+0.0001);
			
// 		_estimatedChannelMatrices[iObservationToBeProcessed] = _channelEstimator->nextMatrix(observations.col(iObservationToBeProcessed),softEstimations,noiseVariances[iObservationToBeProcessed]);
		
		// ...and updated with the last estimated channel matrix
		ARmatricesBuffer.erase(ARmatricesBuffer.begin());
		ARmatricesBuffer.push_back(_estimatedChannelMatrices[iObservationToBeProcessed]);

    } // for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

std::vector<MatrixXd> LinearFilterBasedAlgorithm::getChannelMatricesToStackForSmoothing(std::vector<MatrixXd> ARmatricesBuffer) const
{
	std::vector<double>::const_iterator iterARcoeffs;
	std::vector<MatrixXd>::reverse_iterator iterMatrices;
	
	std::vector<MatrixXd> res(_d+1);
	
	for(uint iSmoothing=0;iSmoothing<(_d+1);iSmoothing++)
	{
		// new matrix to be stacked is initialized to zero
		res[iSmoothing] = MatrixXd::Zero(_nOutputs,_nInputs*(_d+1));
		
		for(iterARcoeffs = _ARcoefficients.begin(),iterMatrices = ARmatricesBuffer.rbegin();iterARcoeffs!=_ARcoefficients.end();iterARcoeffs++,iterMatrices++)
			res[iSmoothing] +=  (*iterARcoeffs)*(*iterMatrices);
		
		ARmatricesBuffer.erase(ARmatricesBuffer.begin());
		ARmatricesBuffer.push_back(res[iSmoothing]);
	}
	
	return res;
}

MatrixXd LinearFilterBasedAlgorithm::getDetectedSymbolVectors()
{
    return _detectedSymbolVectors.block(0,_preamble.cols(),_nInputs,_iLastSymbolVectorToBeDetected-_preamble.cols());
}

vector<MatrixXd> LinearFilterBasedAlgorithm::getEstimatedChannelMatrices()
{
    vector<MatrixXd> channelMatrices;
    channelMatrices.reserve(_iLastSymbolVectorToBeDetected-_preamble.cols());

    for(uint i=_preamble.cols();i<_iLastSymbolVectorToBeDetected;i++)
        channelMatrices.push_back(_estimatedChannelMatrices[i]);

    return channelMatrices;
}
