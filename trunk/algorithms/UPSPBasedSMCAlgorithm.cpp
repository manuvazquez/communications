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
#include "UPSPBasedSMCAlgorithm.h"

// #define DEBUG

UPSPBasedSMCAlgorithm::UPSPBasedSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance): ChannelOrderEstimatorSMCAlgorithm (name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_particleFilter(new ParticleFilter(nParticles)),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_particlesBestChannelOrders(nParticles)
{
}


UPSPBasedSMCAlgorithm::~UPSPBasedSMCAlgorithm()
{
	delete _particleFilter;
}

void UPSPBasedSMCAlgorithm::InitializeParticles()
{
	vector<ChannelMatrixEstimator *> channelEstimatorsClone(_channelEstimators.size());
	for(uint i=0;i<_candidateOrders.size();i++)
		channelEstimatorsClone[i] = _channelEstimators[i]->Clone();

	// we begin with only one particle
	_particleFilter->AddParticle(new ParticleWithChannelEstimation(1.0,_N,_K+_d,channelEstimatorsClone));
	_particleFilter->GetParticle(0)->SetSymbolVectors(tRange(0,_preamble.cols()-1),_preamble);
}

void UPSPBasedSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	uint nSymbolVectors = (int) pow((double)_alphabet.Length(),(double)_N);
	tRange rMaxChannelOrderMinus1FirstColumns(0,_maxOrder-2),rAll;
	vector<tSymbol> testedVector(_N);
	tVector computedObservations(_L),error(_L);
	vector<double> costs(_candidateOrders.size());
	vector<tVector> computedObservationsVector(_candidateOrders.size());
	int iCandidate;
	double normConst;

	typedef struct{
		int fromParticle;
		tMatrix symbolVectorsMatrix;
        int iBestChannelOrder;
        vector<tVector> computedObservationsVector;
		double weight;
	}tParticleCandidate;

	tParticleCandidate *particleCandidates = new tParticleCandidate[_particleFilter->Capacity()*nSymbolVectors];

    double *unnormalizedChannelOrderAPPs = new double[_candidateOrders.size()];
//     double normConst;

    // "symbolVectorsMatrix" will contain all the symbols involved in the current observation
    tMatrix symbolVectorsMatrix(_N,_maxOrder);
    tVector symbolsVector(_NmaxOrder);

    int lastSymbolVectorStart = _NmaxOrder - _N;

	tRange rMaxChannelOrderMinus1PrecedentColumns(_startDetectionTime-_maxOrder+1,_startDetectionTime-1);

	// at first, there is only one particle
	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K+_d;iObservationToBeProcessed++)
	{
		// it keeps track of the place where a new tParticleCandidate will be stored within the array
		iCandidate = 0;

		normConst = 0.0;

		// the candidates from all the particles are generated
		for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
		{
			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			symbolVectorsMatrix(rAll,rMaxChannelOrderMinus1FirstColumns).inject(processedParticle->GetSymbolVectors(rMaxChannelOrderMinus1PrecedentColumns));
			symbolsVector = Util::ToVector(symbolVectorsMatrix,columnwise);

			for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
			{
				// the corresponding testing vector is generated from the index
				_alphabet.IntToSymbolsArray(iTestedVector,testedVector);

				// current tested vector is copied in the m-th position
				for(int k=0;k<_N;k++)
					symbolVectorsMatrix(k,_maxOrder-1) = symbolsVector(lastSymbolVectorStart+k) = testedVector[k];


				for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
				{
					int m = _candidateOrders[iChannelOrder];

					tMatrix estimatedChannelMatrix = processedParticle->GetChannelMatrixEstimator(iChannelOrder)->LastEstimatedChannelMatrix();
					estimatedChannelMatrix *= _ARcoefficient;

					// computedObservations = estimatedChannelMatrix * symbolVectorsMatrix(:)
					Blas_Mat_Vec_Mult(estimatedChannelMatrix,symbolsVector(tRange(_NmaxOrder-m*_N,_NmaxOrder-1)),computedObservations);

					// error = observations.col(iObservationToBeProcessed) - computedObservations
					Util::Add(tVector(observations.col(iObservationToBeProcessed)),tVector(computedObservations),error,1.0,-1.0);

					costs[iChannelOrder] = Blas_Dot_Prod(error,error);
					computedObservationsVector[iChannelOrder] = computedObservations;
				}

				int iLeastCost;
				Util::Min(costs,iLeastCost);

				particleCandidates[iCandidate].fromParticle = iParticle;
				particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;
                particleCandidates[iCandidate].iBestChannelOrder = iLeastCost;
                particleCandidates[iCandidate].computedObservationsVector = computedObservationsVector;
				particleCandidates[iCandidate].weight = processedParticle->GetWeight()*StatUtil::NormalPdf(observations.col(iObservationToBeProcessed),computedObservationsVector[iLeastCost],noiseVariances[iObservationToBeProcessed]);
				normConst += particleCandidates[iCandidate].weight;

				iCandidate++;
			} // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

		} // for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)

		// a vector of size the number of generated candidates is declared...
		tVector weights(iCandidate);

		// ...to store their weights
		for(int i=0;i<iCandidate;i++)
			weights(i) = particleCandidates[i].weight/normConst;

		// the candidates that are going to give particles are selected
		vector<int> indexesSelectedCandidates = _resamplingAlgorithm->ObtainIndexes(_particleFilter->Capacity(),weights);

		// every survivor candidate is associated with an old particle
		vector<int> indexesParticles(indexesSelectedCandidates.size());
		for(uint i=0;i<indexesSelectedCandidates.size();i++)
        {
			indexesParticles[i] = particleCandidates[indexesSelectedCandidates[i]].fromParticle;
            _particlesBestChannelOrders[i] = particleCandidates[indexesSelectedCandidates[i]].iBestChannelOrder;
        }

		// the chosen particles are kept without modification (yet)
		_particleFilter->KeepParticles(indexesParticles);

		// every surviving particle is modified according to what it says its corresponding candidate
		for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)
		{
			ParticleWithChannelEstimation *processedParticle = _particleFilter->GetParticle(iParticle);

			// sampled symbols are copied into the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.col(_maxOrder-1));

			for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				// channel matrix is estimated by means of the particle channel estimator
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed,processedParticle->GetChannelMatrixEstimator(iChannelOrder)->NextMatrix(observations.col(iObservationToBeProcessed),particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix(rAll,tRange(_maxOrder-_candidateOrders[iChannelOrder],_maxOrder-1)),noiseVariances[iObservationToBeProcessed]));
			}

			processedParticle->SetWeight(particleCandidates[indexesSelectedCandidates[iParticle]].weight);

		} // for(int iParticle=0;iParticle<_particleFilter->Nparticles();iParticle++)

		_particleFilter->NormalizeWeights();

        int iBestParticle = _particleFilter->iBestParticle();

        uint iChannelOrder;
        normConst = 0.0;

        for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        {
            unnormalizedChannelOrderAPPs[iChannelOrder] = StatUtil::NormalPdf(observations.col(iObservationToBeProcessed),particleCandidates[indexesSelectedCandidates[iBestParticle]].computedObservationsVector[iChannelOrder],noiseVariances[iObservationToBeProcessed]);
            normConst += unnormalizedChannelOrderAPPs[iChannelOrder];
        }

        if(normConst == 0.0)
            throw RuntimeException("UPSPBasedSMCAlgorithm::Process: normalization constant is zero.");

        for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
            _channelOrderAPPs(iChannelOrder,iObservationToBeProcessed) = unnormalizedChannelOrderAPPs[iChannelOrder]/normConst;

		rMaxChannelOrderMinus1PrecedentColumns = rMaxChannelOrderMinus1PrecedentColumns + 1;

	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K+_d;iObservationToBeProcessed++)
	delete[] particleCandidates;
    delete[] unnormalizedChannelOrderAPPs;
}

int UPSPBasedSMCAlgorithm::BestChannelOrderIndex(int iBestParticle)
{
    return _particlesBestChannelOrders[iBestParticle];
}

vector<vector<tMatrix> > UPSPBasedSMCAlgorithm::ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
    tMatrix sequenceToProcess = Util::Append(_preamble,trainingSequence);

    if(observations.cols() < (_iFirstObservation+trainingSequence.cols()))
        throw RuntimeException("UPSPBasedSMCAlgorithm::ProcessTrainingSequence: Insufficient number of observations.");

	// channel estimation for the training sequence is needed in order to compute the channel order APP
	vector<vector<tMatrix> > estimatedMatrices = UnknownChannelOrderAlgorithm::ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

	tVector computedObservations(_L),unnormalizedChannelOrderAPPs(_candidateOrders.size());
	tRange rAll;

	for(int i=_preamble.cols();i<sequenceToProcess.cols();i++)
	{
		double normConst = 0.0;

		for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
		{
			// computedObservations = estimatedMatrices[iOrder][i-_preamble.cols()] * Util::ToVector(sequenceToProcess(rAll,tRange(i-_candidateOrders[iOrder]+1,i)),columnwise)
			Blas_Mat_Vec_Mult(estimatedMatrices[iOrder][i-_preamble.cols()],Util::ToVector(sequenceToProcess(rAll,tRange(i-_candidateOrders[iOrder]+1,i)),columnwise),computedObservations);

			unnormalizedChannelOrderAPPs(iOrder) = StatUtil::NormalPdf(observations.col(i),computedObservations,noiseVariances[i]);
			normConst += unnormalizedChannelOrderAPPs(iOrder);
		}

		for(uint iOrder=0;iOrder<_candidateOrders.size();iOrder++)
			_channelOrderAPPs(iOrder,i) = unnormalizedChannelOrderAPPs(iOrder)/normConst;
	}

	return estimatedMatrices;
}
