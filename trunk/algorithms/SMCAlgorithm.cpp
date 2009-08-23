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
#include "SMCAlgorithm.h"

// #define DEBUG

SMCAlgorithm::SMCAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m, ChannelMatrixEstimator *channelEstimator, tMatrix preamble,int smoothingLag,int nParticles,ResamplingAlgorithm *resamplingAlgorithm, const tMatrix &channelMatrixMean, const tMatrix &channelMatrixVariances): KnownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected,m, channelEstimator, preamble),
// _variables initialization
_particleFilter(new ParticleFilter(nParticles)),
_particleFilterNeedToBeDeleted(true),_resamplingAlgorithm(resamplingAlgorithm),_d(smoothingLag),_estimatorIndex(0),
_channelMean(Util::toVector(Util::lapack2eigen(channelMatrixMean),rowwise)),_channelCovariance(Util::toVector(Util::lapack2eigen(channelMatrixVariances),rowwise).asDiagonal()),_randomParticlesInitilization(false)
{
    if(channelMatrixMean.rows()!=Nr || channelMatrixMean.cols()!=(N*m))
    {
        cout << "channelMatrixMean.rows() = " << channelMatrixMean.rows() << " channelMatrixMean.cols() = " << channelMatrixMean.cols() << endl;
        throw RuntimeException("SMCAlgorithm::SMCAlgorithm: channel matrix mean dimensions are wrong.");
    }

    if(channelMatrixVariances.rows()!=Nr || channelMatrixVariances.cols()!=(N*m))
        throw RuntimeException("SMCAlgorithm::SMCAlgorithm: channel matrix variances dimensions are wrong.");

    // at first, we assume that all observations from the preamble need to be processed
    _startDetectionTime = _preamble.cols();
}

// constructor that receives an already functional particle filter
SMCAlgorithm::SMCAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m, tMatrix preamble,int smoothingLag,ParticleFilter *particleFilter,ResamplingAlgorithm *resamplingAlgorithm): KnownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected,m, preamble),
// _variables initialization
_particleFilter(particleFilter),_particleFilterNeedToBeDeleted(false),_resamplingAlgorithm(resamplingAlgorithm),_d(smoothingLag),_estimatorIndex(0)
{
    // at first, we assume that all observations from the preamble need to be processed
    _startDetectionTime = _preamble.cols();
}

SMCAlgorithm::~SMCAlgorithm()
{
    if(_particleFilterNeedToBeDeleted)
        delete _particleFilter;
}

void SMCAlgorithm::setEstimatorIndex(int n)
{
    if(_particleFilter==NULL)
        throw RuntimeException("SMCAlgorithm::SetEstimatorIndex: the particle filter is not set.");

    if(n>=dynamic_cast<WithChannelEstimationParticleAddon *>(_particleFilter->getParticle(0))->nChannelMatrixEstimators())
        throw RuntimeException("SMCAlgorithm::SetEstimatorIndex: index is out of range.");

    _estimatorIndex = n;
}

void SMCAlgorithm::initializeParticles()
{
    ChannelMatrixEstimator *channelMatrixEstimatorClone;

    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
    {
        channelMatrixEstimatorClone = _channelEstimator->clone();
        if(_randomParticlesInitilization)
            channelMatrixEstimatorClone->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::randnMatrix(_channelMean,_channelCovariance),rowwise,_Nr));
        _particleFilter->addParticle(new ParticleWithChannelEstimation(1.0/(double)_particleFilter->capacity(),_nInputs,_iLastSymbolVectorToBeDetected,channelMatrixEstimatorClone));

        // if there is preamble...
        if(_preamble.cols()!=0)
            _particleFilter->getParticle(iParticle)->setSymbolVectors(0,_preamble.cols(),Util::lapack2eigen(_preamble));
    }
}

void SMCAlgorithm::run(MatrixXd observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();

    if(nObservations<(_startDetectionTime+1+_d))
        throw RuntimeException("SMCAlgorithm::run: not enough observations.");

    initializeParticles();

    process(observations,noiseVariances);
}

void SMCAlgorithm::runFrom(int n,MatrixXd observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();
    _startDetectionTime = n;

    if(nObservations<(_startDetectionTime+1+_d))
        throw RuntimeException("SMCAlgorithm::runFrom: Not enough observations.");

    process(observations,noiseVariances);
}

void SMCAlgorithm::run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence)
{
    if(observations.rows()!=_nOutputs || trainingSequence.rows()!=_nInputs)
        throw RuntimeException("SMCAlgorithm::Run: Observations matrix or training sequence dimensions are wrong.");

    int iParticle,j;
    
    MatrixXd preambleTrainingSequence(trainingSequence.rows(),_preamble.cols()+trainingSequence.cols());
    preambleTrainingSequence << Util::lapack2eigen(_preamble),trainingSequence;    
    
    beforeInitializingParticles(observations,trainingSequence);

    initializeParticles();

    for(iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)
    {
        WithChannelEstimationParticleAddon *processedParticle = dynamic_cast<WithChannelEstimationParticleAddon *>(_particleFilter->getParticle(iParticle));

        vector<MatrixXd> trainingSequenceChannelMatrices = processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatricesFromObservationsSequence(observations,noiseVariances,preambleTrainingSequence,_preamble.cols(),preambleTrainingSequence.cols());

        //the channel estimation given by the training sequence is copied into each particle...
        for(j=_preamble.cols();j<preambleTrainingSequence.cols();j++)
            processedParticle->setChannelMatrix(_estimatorIndex,j,trainingSequenceChannelMatrices[j-_preamble.cols()]);

        //... the symbols are considered detected...
        _particleFilter->getParticle(iParticle)->setSymbolVectors(_preamble.cols(),preambleTrainingSequence.cols(),trainingSequence);
    }

    // the process method must start in
    _startDetectionTime = preambleTrainingSequence.cols();

    process(observations,noiseVariances);
}

MatrixXd SMCAlgorithm::getDetectedSymbolVectors_eigen()
{
    return _particleFilter->getBestParticle()->getSymbolVectors_eigen(_preamble.cols(),_iLastSymbolVectorToBeDetected-1);
}

vector<MatrixXd> SMCAlgorithm::getEstimatedChannelMatrices_eigen()
{
    vector<MatrixXd> channelMatrices;
    channelMatrices.reserve(_iLastSymbolVectorToBeDetected-_preamble.cols());

    // best particle is chosen
    int iBestParticle = _particleFilter->iBestParticle();

    for(int i=_preamble.cols();i<_iLastSymbolVectorToBeDetected;i++)
        channelMatrices.push_back(dynamic_cast<WithChannelEstimationParticleAddon *>(_particleFilter->getParticle(iBestParticle))->getChannelMatrix_eigen(_estimatorIndex,i));

    return channelMatrices;
}

double SMCAlgorithm::smoothedLikelihood(const vector<MatrixXd> &channelMatrices,const MatrixXd &involvedSymbolVectors,int iObservationToBeProcessed,const MatrixXd &observations,const vector<double> &noiseVariances)
{
    double likelihoodsProd = 1.0;

    for(int iSmoothing=0;iSmoothing<=_d;iSmoothing++)
    {
        VectorXd stackedSymbolVector = Util::toVector(involvedSymbolVectors.block(0,iSmoothing,_nInputs,_channelOrder),columnwise);

        likelihoodsProd *= StatUtil::normalPdf(observations.col(iObservationToBeProcessed+iSmoothing),channelMatrices[iSmoothing]*stackedSymbolVector,noiseVariances[iObservationToBeProcessed+iSmoothing]);
    }
    return likelihoodsProd;
}
