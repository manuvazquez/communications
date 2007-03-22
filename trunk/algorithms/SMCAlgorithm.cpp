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

// #define DEBUG10

SMCAlgorithm::SMCAlgorithm(string name, Alphabet alphabet,int L,int N, int K,int m, ChannelMatrixEstimator *channelEstimator, tMatrix preamble,int smoothingLag,int nParticles,ResamplingAlgorithm *resamplingAlgorithm, const tMatrix &channelMatrixMean, const tMatrix &channelMatrixVariances): KnownChannelOrderAlgorithm(name, alphabet, L, N, K,m, channelEstimator, preamble),
// _variables initialization
_particleFilter(new ParticleFilter(nParticles)),_particleFilterNeedToBeDeleted(true),_resamplingAlgorithm(resamplingAlgorithm),_d(smoothingLag),_allSymbolsRows(0,_N-1),_estimatorIndex(0),_channelMatrixMean(channelMatrixMean),_channelMatrixVariances(channelMatrixVariances)
{
	if(channelMatrixMean.rows()!=L || channelMatrixMean.cols()!=(N*m))
		throw RuntimeException("SMCAlgorithm::SMCAlgorithm: channel matrix mean dimensions are wrong.");

	if(channelMatrixVariances.rows()!=L || channelMatrixVariances.cols()!=(N*m))
		throw RuntimeException("SMCAlgorithm::SMCAlgorithm: channel matrix variances dimensions are wrong.");

    // at first, we assume that all observations from the preamble need to be processed
    _startDetectionTime = _preamble.cols();
}

SMCAlgorithm::SMCAlgorithm(string name, Alphabet alphabet,int L,int N, int K,int m, tMatrix preamble,int smoothingLag,ParticleFilter *particleFilter,ResamplingAlgorithm *resamplingAlgorithm): KnownChannelOrderAlgorithm(name, alphabet, L, N, K,m, preamble),
// _variables initialization
_particleFilter(particleFilter),_particleFilterNeedToBeDeleted(false),_resamplingAlgorithm(resamplingAlgorithm),_d(smoothingLag),_allSymbolsRows(0,_N-1),_estimatorIndex(0)
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
    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
    {
        _particleFilter->SetParticle(new ParticleWithChannelEstimation(1.0/(double)_particleFilter->Nparticles(),_N,_K,_channelEstimator->Clone()),iParticle);

        _particleFilter->GetParticle(iParticle)->SetSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);

        #ifdef DEBUG10
            cout << "lo que acabo de meter" << endl << _particleFilter->GetParticle(iParticle)->GetSymbolVectors(tRange(0,_preamble.cols()-1));
        #endif
    }
}

void SMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();

    if(nObservations<(_startDetectionTime+1+_d))
        throw RuntimeException("SMCAlgorithm::Run: Not enough observations.");

    this->InitializeParticles();

	tVector channelMean = Util::ToVector(_channelMatrixMean,rowwise);
	tMatrix channelCovariance = LaGenMatDouble::from_diag(Util::ToVector(_channelMatrixVariances,rowwise));

	#ifdef DEBUG
		cout << "La media es" << endl << channelMean;
		cout << "La covarianza es" << endl << channelCovariance;
	#endif

	int nDimensions = _L*_Nm;
	int nPoinsByDimension = pow((double)_particleFilter->Nparticles(),1.0/(double)nDimensions);
	int remainingParticles = _particleFilter->Nparticles() - (int)pow((double)nPoinsByDimension,(double)nDimensions);

	#ifdef DEBUG
		cout << "Sobran " << remainingParticles << " partículas" << endl;
	#endif

	vector<vector<double> > alphabets(nDimensions,vector<double>(nPoinsByDimension));
	for(int i=0;i<nDimensions;i++)
	{
		double sigma = sqrt(channelCovariance(i,i));
		double step = 4.0*sigma/(double)(nPoinsByDimension+1);
		double leftBound = channelMean(i) -  2.0*sigma;
		for(int j=0;j<nPoinsByDimension;j++)
		{
			leftBound += step;
			alphabets[i][j] = leftBound;
		}
	}

	#ifdef DEBUG
		cout << "Los alfabetos son" << endl;
		for(int i=0;i<nDimensions;i++)
		{
			Util::Print(alphabets[i]);
		}
	#endif


	vector<double> channelSample(nDimensions);
	for(int i=0;i<nDimensions;i++)
		channelSample[i] = alphabets[i][0];

	#ifdef DEBUG
		cout << "nPoinsByDimension = " << nPoinsByDimension << endl;
		exit(0);
	#endif

	// the initial estimation of the particles channel matrix estimators is set
    for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
    {
		ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

		// a sample of the a priori channel distribution is drawn based on the mean and the covariance
		processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->SetFirstEstimatedChannelMatrix(Util::ToMatrix( channelSample,rowwise,_L,_Nm));
		Util::NextVector(channelSample,alphabets);

// 		processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->SetFirstEstimatedChannelMatrix(Util::ToMatrix( StatUtil::RandMatrix(channelMean,channelCovariance),rowwise,_L));
    }

    this->Process(observations,noiseVariances);
}

void SMCAlgorithm::RunFrom(int n,tMatrix observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();
	_startDetectionTime = n;

    if(nObservations<(_startDetectionTime+1+_d))
        throw RuntimeException("SMCAlgorithm::RunFrom: Not enough observations.");

    this->Process(observations,noiseVariances);
}

void SMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
    if(observations.rows()!=_L || trainingSequence.rows()!=_N)
        throw RuntimeException("SMCAlgorithm::Run: Observations matrix or training sequence dimensions are wrong.");

    int iParticle,j;
    int preamblePlusTrainingSequenceLength = _preamble.cols() + trainingSequence.cols();

    tRange rTrainingSequence(_preamble.cols(),preamblePlusTrainingSequenceLength-1);

    vector<tMatrix> trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

    this->InitializeParticles();

    for(iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
    {
		ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

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

    this->Process(observations,noiseVariances);
}

tMatrix SMCAlgorithm::GetDetectedSymbolVectors()
{
    // best particle is chosen
    int iBestParticle;
    Util::Max(_particleFilter->GetWeightsVector(),iBestParticle);

    return ((_particleFilter->GetParticle(iBestParticle))->GetAllSymbolVectors())(_allSymbolsRows,tRange(_preamble.cols(),_K-1));
}

vector<tMatrix> SMCAlgorithm::GetEstimatedChannelMatrices()
{
    vector<tMatrix> channelMatrices;
    channelMatrices.reserve(_K-_preamble.cols());

    // best particle is chosen
    int iBestParticle;
    Util::Max(_particleFilter->GetWeightsVector(),iBestParticle);

    for(int i=_preamble.cols();i<_K;i++)
        channelMatrices.push_back((_particleFilter->GetParticle(iBestParticle))->GetChannelMatrix(_estimatorIndex,i));

    return channelMatrices;
}
