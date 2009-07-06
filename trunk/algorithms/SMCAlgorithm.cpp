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
_particleFilterNeedToBeDeleted(true),_resamplingAlgorithm(resamplingAlgorithm),_d(smoothingLag),_allSymbolsRows(0,_nInputs-1),_estimatorIndex(0),
_channelMatrixMean(channelMatrixMean),_channelMatrixVariances(channelMatrixVariances),_randomParticlesInitilization(false)
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
_particleFilter(particleFilter),_particleFilterNeedToBeDeleted(false),_resamplingAlgorithm(resamplingAlgorithm),_d(smoothingLag),_allSymbolsRows(0,_nInputs-1),_estimatorIndex(0)
{
    // at first, we assume that all observations from the preamble need to be processed
    _startDetectionTime = _preamble.cols();
}

SMCAlgorithm::~SMCAlgorithm()
{
    if(_particleFilterNeedToBeDeleted)
        delete _particleFilter;
}

void SMCAlgorithm::SetEstimatorIndex(int n)
{
    if(_particleFilter==NULL)
        throw RuntimeException("SMCAlgorithm::SetEstimatorIndex: the particle filter is not set.");

    if(n>=_particleFilter->GetParticle(0)->NchannelMatrixEstimators())
        throw RuntimeException("SMCAlgorithm::SetEstimatorIndex: index is out of range.");

    _estimatorIndex = n;
}

void SMCAlgorithm::InitializeParticles()
{
    ChannelMatrixEstimator *channelMatrixEstimatorClone;
    tVector channelMean = Util::toVector(_channelMatrixMean,rowwise);
    tMatrix channelCovariance = LaGenMatDouble::from_diag(Util::toVector(_channelMatrixVariances,rowwise));

    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)
    {
        channelMatrixEstimatorClone = _channelEstimator->Clone();
        if(_randomParticlesInitilization)
            channelMatrixEstimatorClone->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::RandMatrix(channelMean,channelCovariance),rowwise,_nOutputs));
        _particleFilter->AddParticle(new ParticleWithChannelEstimation(1.0/(double)_particleFilter->Capacity(),_nInputs,_iLastSymbolVectorToBeDetected,channelMatrixEstimatorClone));

        // if there is preamble...
        if(_preamble.cols()!=0)
            _particleFilter->GetParticle(iParticle)->SetSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);
    }
}

void SMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();

    if(nObservations<(_startDetectionTime+1+_d))
        throw RuntimeException("SMCAlgorithm::Run: not enough observations.");

    InitializeParticles();

    Process(observations,noiseVariances);
}

void SMCAlgorithm::RunFrom(int n,tMatrix observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();
    _startDetectionTime = n;

    if(nObservations<(_startDetectionTime+1+_d))
        throw RuntimeException("SMCAlgorithm::RunFrom: Not enough observations.");

    Process(observations,noiseVariances);
}

void SMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
    if(observations.rows()!=_nOutputs || trainingSequence.rows()!=_nInputs)
        throw RuntimeException("SMCAlgorithm::Run: Observations matrix or training sequence dimensions are wrong.");

    int iParticle,j;
    tMatrix preamblePlusTrainingSequence = Util::append(_preamble,trainingSequence);
    int preamblePlusTrainingSequenceLength = preamblePlusTrainingSequence.cols();

    tRange rTrainingSequence(_preamble.cols(),preamblePlusTrainingSequenceLength-1);

    BeforeInitializingParticles(observations,trainingSequence);

    InitializeParticles();

//     if(_randomParticlesInitilization)
//         cout << _name << ": particles are initialized randomly..." << endl;

    for(iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
    {
        ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

        vector<tMatrix> trainingSequenceChannelMatrices = processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->nextMatricesFromObservationsSequence(observations,noiseVariances,preamblePlusTrainingSequence,_preamble.cols(),preamblePlusTrainingSequenceLength);

        //the channel estimation given by the training sequence is copied into each particle...
        for(j=_preamble.cols();j<preamblePlusTrainingSequenceLength;j++)
        {
            processedParticle->SetChannelMatrix(_estimatorIndex,j,trainingSequenceChannelMatrices[j-_preamble.cols()]);
        }

        //... the symbols are considered detected...
        processedParticle->SetSymbolVectors(rTrainingSequence,trainingSequence);
    }

    // the Process method must start in
    _startDetectionTime = preamblePlusTrainingSequenceLength;

    Process(observations,noiseVariances);
}

tMatrix SMCAlgorithm::getDetectedSymbolVectors()
{
    return (_particleFilter->GetBestParticle()->GetAllSymbolVectors())(_allSymbolsRows,tRange(_preamble.cols(),_iLastSymbolVectorToBeDetected-1));
}

vector<tMatrix> SMCAlgorithm::GetEstimatedChannelMatrices()
{
    vector<tMatrix> channelMatrices;
    channelMatrices.reserve(_iLastSymbolVectorToBeDetected-_preamble.cols());

    // best particle is chosen
    int iBestParticle = _particleFilter->iBestParticle();

    for(int i=_preamble.cols();i<_iLastSymbolVectorToBeDetected;i++)
        channelMatrices.push_back((_particleFilter->GetParticle(iBestParticle))->GetChannelMatrix(_estimatorIndex,i));

    return channelMatrices;
}

double SMCAlgorithm::smoothedLikelihood(const vector<tMatrix> &channelMatrices,const tMatrix &involvedSymbolVectors,ParticleWithChannelEstimation *particle,int iObservationToBeProcessed,const tMatrix &observations,const vector<double> &noiseVariances)
{
    double likelihoodsProd = 1.0;

    tVector predictedNoiselessObservation(_nOutputs);
    tRange rAllSymbolRows(0,_nInputs-1);

    for(int iSmoothing=0;iSmoothing<=_d;iSmoothing++)
    {
        tRange rSymbolVectors(iSmoothing,iSmoothing+_channelOrder-1);
        tVector stackedSymbolVector = Util::toVector(involvedSymbolVectors(rAllSymbolRows,rSymbolVectors),columnwise);

        // predictedNoiselessObservation = matricesToStack[iSmoothing] * stackedSymbolVector
        Blas_Mat_Vec_Mult(channelMatrices[iSmoothing],stackedSymbolVector,predictedNoiselessObservation);

        likelihoodsProd *= StatUtil::NormalPdf(observations.col(iObservationToBeProcessed+iSmoothing),predictedNoiselessObservation,noiseVariances[iObservationToBeProcessed+iSmoothing]);
    }
    return likelihoodsProd;
}
