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
#include "TriangularizationBasedSMCAlgorithm.h"

// #define DEBUG4

TriangularizationBasedSMCAlgorithm::TriangularizationBasedSMCAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm, const MatrixXd& channelMatrixMean, const MatrixXd& channelMatrixVariances,double ARcoefficient,double ARprocessVariance): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances),_ARcoefficient(ARcoefficient),_ARprocessVariance(ARprocessVariance)
{
//     _randomParticlesInitilization = true;
}

void TriangularizationBasedSMCAlgorithm::process(const MatrixXd& observations, vector<double> noiseVariances)
{
    uint iParticle,iSmoothing,iAlphabet,iSampled;
    double proposal,observationWithouNoise,sumProb,likelihoodsProd;
    vector<MatrixXd> matricesToStack(_d+1,MatrixXd(_nOutputs,_nInputsXchannelOrder));
    MatrixXd observationsCovariance = MatrixXd::Zero(_nOutputs*(_d+1),_nOutputs*(_d+1));
    MatrixXd involvedSymbolVectors = MatrixXd::Zero(_nInputs,_channelOrder+_d);
    uint NmMinus1 = _nInputs*(_channelOrder-1);
    VectorXd symbolProbabilities(_alphabet.length());

    for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
        // the stacked observations vector
        VectorXd stackedObservations = Util::toVector(observations.block(0,iObservationToBeProcessed,_nOutputs,_d+1),columnwise);

        for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
        {
            ParticleWithChannelEstimation *processedParticle = dynamic_cast <ParticleWithChannelEstimation *> (_particleFilter->getParticle(iParticle));

            // the already detected symbol vectors are stored in "involvedSymbolVectors"
            involvedSymbolVectors.block(0,0,_nInputs,_channelOrder-1) = processedParticle->getSymbolVectors(iObservationToBeProcessed-_channelOrder+1,iObservationToBeProcessed-1);

            // predicted channel matrices are stored in a vector in order to stack them
            // (first one is obtained via the Kalman Filter)
            matricesToStack[0] = (dynamic_cast<KalmanEstimator *> (processedParticle->getChannelMatrixEstimator(_estimatorIndex)))->samplePredicted();

            for(iSmoothing=1;iSmoothing<=_d;iSmoothing++)
                matricesToStack[iSmoothing] = _ARcoefficient*matricesToStack[iSmoothing-1]+StatUtil::randnMatrix(_nOutputs,_nInputsXchannelOrder,0.0,_ARprocessVariance);

            // matrices are stacked to give
            MatrixXd stackedChannelMatrix = channelMatrices2stackedChannelMatrix(matricesToStack);

            MatrixXd stackedChannelMatrixMinus = stackedChannelMatrix.block(0,(_channelOrder-1)*_nInputs,stackedChannelMatrix.rows(),stackedChannelMatrix.cols()-(_channelOrder-1)*_nInputs);

            VectorXd stackedObservationsMinus = substractKnownSymbolsContribution(matricesToStack,_channelOrder,_d,stackedObservations,involvedSymbolVectors.block(0,0,_nInputs,_channelOrder-1));
            
            // we want to start sampling the present symbol vector, not the future ones
            MatrixXd stackedChannelMatrixMinusFlipped = Util::flipLR(stackedChannelMatrixMinus);

            Eigen::LLT<MatrixXd> llt(stackedChannelMatrixMinusFlipped.transpose()*stackedChannelMatrixMinusFlipped);

			MatrixXd U = llt.matrixL().transpose();

			MatrixXd invLstackedChannelMatrixMinusTrans = llt.matrixL().solve(stackedChannelMatrixMinusFlipped.transpose());
            VectorXd transformedStackedObservationsMinus = invLstackedChannelMatrixMinusTrans*stackedObservationsMinus;

            // the covariance of the transformed observations is computed...

            // ...starting by the covariance of the normal observations
            for(iSmoothing=0;iSmoothing<_d+1;iSmoothing++)
                for(uint i=0;i<_nOutputs;i++)
                    observationsCovariance(iSmoothing*_nOutputs+i,iSmoothing*_nOutputs+i) = noiseVariances[iObservationToBeProcessed+iSmoothing];

            MatrixXd transformedStackedObservationsCovariance = invLstackedChannelMatrixMinusTrans*observationsCovariance*invLstackedChannelMatrixMinusTrans.transpose();

            // the evaluated proposal function (necessary for computing the weights) is initialized
            proposal = 1.0;

            for(int iSampledSymbol=_nInputs*(_d+1)-1,iWithinMatrix=NmMinus1;iSampledSymbol>=0;iSampledSymbol--,iWithinMatrix++)
            {
                observationWithouNoise = 0.0;
                int jU,iS;
                for(jU=_nInputs*(_d+1)-1,iS=NmMinus1;jU>iSampledSymbol;jU--,iS++)
                    observationWithouNoise += U(iSampledSymbol,jU)*involvedSymbolVectors(iS % _nInputs,iS / _nInputs);

                sumProb = 0.0;

                // the probability for each posible symbol alphabet is computed
                for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                {
                    symbolProbabilities(iAlphabet) = StatUtil::normalPdf(transformedStackedObservationsMinus(iSampledSymbol),observationWithouNoise+U(iSampledSymbol,jU)*_alphabet[iAlphabet],transformedStackedObservationsCovariance(iSampledSymbol,iSampledSymbol));

                    // for normalization purposes
                    sumProb += symbolProbabilities(iAlphabet);
                }

                if(sumProb!=0)
                    for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                        symbolProbabilities(iAlphabet) /= sumProb;
                else
                {                        
                    cout << "TriangularizationBasedSMCAlgorithm::process: the sum of the probabilities is null." << endl;
                    cout <<  __FILE__  << "(line " << __LINE__ << ") :" << endl;
                    for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                        symbolProbabilities(iAlphabet) = 1.0/double(_alphabet.length());
                }

                iSampled = StatUtil::discrete_rnd(symbolProbabilities);
                involvedSymbolVectors(iWithinMatrix % _nInputs,iWithinMatrix / _nInputs) = _alphabet[iSampled];
                proposal *= symbolProbabilities(iSampled);

            } // for(int iSampledSymbol=_nInputs*(_d+1)-1,iWithinMatrix=NmMinus1;iSampledSymbol>=0;iSampledSymbol--,iWithinMatrix++)

            processedParticle->setSymbolVector(iObservationToBeProcessed,involvedSymbolVectors.col(_channelOrder-1));

            likelihoodsProd = smoothedLikelihood(matricesToStack,involvedSymbolVectors,iObservationToBeProcessed,observations,noiseVariances);

            // the weight is updated
            processedParticle->setWeight((likelihoodsProd/proposal)*processedParticle->getWeight());

            // and the estimation of the channel matrix
            processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,
                                                processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),
                                                involvedSymbolVectors.block(0,0,_nInputs,_channelOrder),noiseVariances[iObservationToBeProcessed]));            
        } // for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)

        _particleFilter->normalizeWeights();

        // if it's not the last time instant
        if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
            _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);

    } // for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}
