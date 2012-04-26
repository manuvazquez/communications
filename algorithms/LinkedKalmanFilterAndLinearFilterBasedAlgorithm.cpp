/*
    Copyright 2012 Manu <manuavazquez@gmail.com>

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/


#include "LinkedKalmanFilterAndLinearFilterBasedAlgorithm.h"

LinkedKalmanFilterAndLinearFilterBasedAlgorithm::LinkedKalmanFilterAndLinearFilterBasedAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator* channelEstimator, MatrixXd preamble, uint smoothingLag, LinearDetector *linearDetector,  double ARcoefficient, bool substractContributionFromKnownSymbols)
:LinearFilterBasedAlgorithm(name,alphabet,L,Nr,N,iLastSymbolVectorToBeDetected,m,channelEstimator,preamble,smoothingLag,linearDetector,ARcoefficient,substractContributionFromKnownSymbols)
{
	(dynamic_cast<LinearFilterAwareNoiseVarianceAdjustingKalmanEstimatorDecorator *> (_channelEstimator))->setLinearDetector(_linearDetector);
}

// void LinkedKalmanFilterAndLinearFilterBasedAlgorithm::process(const MatrixXd &observations,vector<double> noiseVariances, MatrixXd trainingSequence)
// {
// 	assert(observations.rows()==_nOutputs); 
// 	
//     uint nObservations = observations.cols();
// 
//     uint startDetectionTime = _preamble.cols() + trainingSequence.cols();
// 
// 	assert(nObservations>=startDetectionTime+1+_d);
// 
//     if(_substractContributionFromKnownSymbols)
//     {
//         if(_linearDetector->channelMatrixcols() != _nInputs*(_d+1))
//             throw RuntimeException("LinearFilterBasedAlgorithm::process: the algorithm is supposed to operate substracting the contribution of the known symbols but this is not compatible with the current linear detector.");
//     }
// 
//     vector<MatrixXd> matricesToStack(_d+1);
//     uint iSmoothing,iRow;
//     MatrixXd stackedNoiseCovariance = MatrixXd::Zero(_nOutputs*(_d+1),_nOutputs*(_d+1));
//     double ARcoefficientPower;
// 
//     for(uint iObservationToBeProcessed=startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
//     {
//         ARcoefficientPower = _ARcoefficient;
//         for(iSmoothing=0;iSmoothing<=_d;iSmoothing++)
//         {
// 			// if there is no training sequence...
// 			if(iObservationToBeProcessed==_preamble.cols())
// 				// the "last estimated channel matrix" given by the channel estimator is used (it should be the matrix that served to initialize the channel estimator)
// 				matricesToStack[iSmoothing] = ARcoefficientPower*_channelEstimator->lastEstimatedChannelMatrix();
// 			else
// 				matricesToStack[iSmoothing] = ARcoefficientPower*_estimatedChannelMatrices[iObservationToBeProcessed-1];
//             ARcoefficientPower *= _ARcoefficient;
//         }
// 
//         MatrixXd stackedChannelMatrix = channelMatrices2stackedChannelMatrix(matricesToStack);
// 
//         VectorXd stackedObservations = Util::toVector(observations.block(0,iObservationToBeProcessed,_nOutputs,_d+1),columnwise);
// 
//         // stacked noise covariance needs to be constructed
//         for(iSmoothing=0;iSmoothing<=_d;iSmoothing++)
//             for(iRow=0;iRow<_nOutputs;iRow++)
//                 stackedNoiseCovariance((iSmoothing)*_nOutputs+iRow,(iSmoothing)*_nOutputs+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];
// 
//         VectorXd softEstimations;
// 
//         if(_substractContributionFromKnownSymbols)
//         {
//             softEstimations =  _linearDetector->detect(
//                 // the last range chooses all the already detected symbol vectors
//                 substractKnownSymbolsContribution(matricesToStack,_channelOrder,_d,stackedObservations,_detectedSymbolVectors.block(0,iObservationToBeProcessed-_channelOrder+1,_nInputs,_channelOrder-1)),
//                 // only a part of the channel matrix is needed. The first range chooses all the stacked observation rows
//                 stackedChannelMatrix.block(0,(_channelOrder-1)*_nInputs,_nOutputs*(_d+1),stackedChannelMatrix.cols()-(_channelOrder-1)*_nInputs),
//                 stackedNoiseCovariance);                
//         } else
//             softEstimations =  _linearDetector->detect(stackedObservations,stackedChannelMatrix,stackedNoiseCovariance);
// 
//         for(iRow=0;iRow<_nInputs;iRow++)
//             _detectedSymbolVectors(iRow,iObservationToBeProcessed) = _alphabet.hardDecision(softEstimations(iRow));
// 
//         _estimatedChannelMatrices[iObservationToBeProcessed] = _channelEstimator->nextMatrix(observations.col(iObservationToBeProcessed),_detectedSymbolVectors.block(0,iObservationToBeProcessed-_channelOrder+1,_nInputs,_channelOrder),noiseVariances[iObservationToBeProcessed]);
// 		
// //         _estimatedChannelMatrices[iObservationToBeProcessed] = _channelEstimator->nextMatrix(observations.col(iObservationToBeProcessed),softEstimations,noiseVariances[iObservationToBeProcessed]);
// 
//     } // for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
// }
