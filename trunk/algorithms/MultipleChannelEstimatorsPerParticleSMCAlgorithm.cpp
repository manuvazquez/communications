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
#include "MultipleChannelEstimatorsPerParticleSMCAlgorithm.h"

// #define DEBUG

MultipleChannelEstimatorsPerParticleSMCAlgorithm::MultipleChannelEstimatorsPerParticleSMCAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation,int smoothingLag,int nParticles,ResamplingAlgorithm *resamplingAlgorithm): UnknownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, preamble, iFirstObservation)
//variables initialization
,_resamplingAlgorithm(resamplingAlgorithm),_d(smoothingLag),_randomParticlesInitilization(false)
{
    // at first, we assume that all symbol vectors from the preamble need to be processed
    _startDetectionTime = _preamble.cols();

    _channelUniqueMean = 0.0;
    _channelUniqueVariance = 1.0;

    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
    {
        _channelMeanVectors.push_back(VectorXd::Constant(_nOutputs*_nInputs*_candidateOrders[iChannelOrder],_channelUniqueMean));
        _channelCovariances.push_back(VectorXd::Constant(_nOutputs*_nInputs*_candidateOrders[iChannelOrder],_channelUniqueMean).asDiagonal());        
    }
}

void MultipleChannelEstimatorsPerParticleSMCAlgorithm::run(MatrixXd observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();

    if(nObservations<_startDetectionTime+_maxOrder)
        throw RuntimeException("MultipleChannelEstimatorsPerParticleSMCAlgorithm::run: not enough observations.");

    this->initializeParticles();

    this->process(observations,noiseVariances);
}

void MultipleChannelEstimatorsPerParticleSMCAlgorithm::run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence)
{
    if(observations.rows()!=_nOutputs || trainingSequence.rows()!=_nInputs)
        throw RuntimeException("MultipleChannelEstimatorsPerParticleSMCAlgorithm::run: observations matrix or training sequence dimensions are wrong.");

    int iParticle;
    uint j;
    uint iChannelOrder;

    // needed for updateParticleChannelOrderEstimators in the USIS algorithm
    vector<vector<MatrixXd> > channelOrderTrainingSequenceChannelMatrices(_candidateOrders.size());

    // to process the training sequence, we need both the preamble and the symbol vectors related to it
    MatrixXd preambleTrainingSequence(trainingSequence.rows(),_preamble.cols()+trainingSequence.cols());
    preambleTrainingSequence << _preamble,trainingSequence;

    tRange rSymbolVectorsTrainingSequece(0,preambleTrainingSequence.cols()-1);

    beforeInitializingParticles(observations,noiseVariances,trainingSequence);

    this->initializeParticles();

    for(iParticle=0;iParticle<getParticleFilterPointer()->nParticles();iParticle++)
    {
        WithChannelEstimationParticleAddon *processedParticle = dynamic_cast <WithChannelEstimationParticleAddon *> (getParticleFilterPointer()->getParticle(iParticle));

        for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        {
            vector<MatrixXd> trainingSequenceChannelMatrices = processedParticle->getChannelMatrixEstimator(iChannelOrder)->nextMatricesFromObservationsSequence(observations,noiseVariances,preambleTrainingSequence,_preamble.cols(),preambleTrainingSequence.cols());
            channelOrderTrainingSequenceChannelMatrices[iChannelOrder] = trainingSequenceChannelMatrices;

            //the channel estimation given by the training sequence is copied into each particle...
            for(j=0;j<trainingSequenceChannelMatrices.size();j++)
                processedParticle->setChannelMatrix(iChannelOrder,_preamble.cols()+j,trainingSequenceChannelMatrices[j]);
        }

        updateParticleChannelOrderEstimators(getParticleFilterPointer()->getParticle(iParticle),observations,channelOrderTrainingSequenceChannelMatrices,noiseVariances,preambleTrainingSequence);

        //... the symbols are considered detected...
        getParticleFilterPointer()->getParticle(iParticle)->setSymbolVectors(0,preambleTrainingSequence.cols(),preambleTrainingSequence);
    }

    // the process method must start in
    _startDetectionTime += trainingSequence.cols();

    this->process(observations,noiseVariances);
}

MatrixXd MultipleChannelEstimatorsPerParticleSMCAlgorithm::getDetectedSymbolVectors_eigen()
{
    return getParticleFilterPointer()->getBestParticle()->getSymbolVectors_eigen(_preamble.cols(),_iLastSymbolVectorToBeDetected-1);
}

vector<MatrixXd> MultipleChannelEstimatorsPerParticleSMCAlgorithm::getEstimatedChannelMatrices_eigen()
{
    vector<MatrixXd> channelMatrices;
    channelMatrices.reserve(_iLastSymbolVectorToBeDetected-_preamble.cols());

    // best particle is chosen
    int iBestParticle = getParticleFilterPointer()->iBestParticle();

    int indexBestChannelOrder = iBestChannelOrder(iBestParticle);

    for(int i=_preamble.cols();i<_iLastSymbolVectorToBeDetected;i++)
        channelMatrices.push_back(dynamic_cast<WithChannelEstimationParticleAddon *>(getParticleFilterPointer()->getParticle(iBestParticle))->getChannelMatrix_eigen(indexBestChannelOrder,i));

    return channelMatrices;
}
