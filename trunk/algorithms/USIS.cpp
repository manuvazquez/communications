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
#include "USIS.h"

// #define DEBUG

USIS::USIS(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators,vector<LinearDetector *> linearDetectors, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,ChannelOrderEstimator * channelOrderEstimator,double ARcoefficient,double samplingVariance,double ARprocessVariance): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_linearDetectors(linearDetectors.size()),_channelOrderEstimator(channelOrderEstimator->Clone()),_particleFilter(nParticles),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_rAllObservationRows(0,_L-1)
#ifdef CHANNELORDERSAPP_SAVING
	,_channelOrderAPPs(LaGenMatDouble::zeros(_candidateOrders.size(),_K)),_rCandidateOrders(0,_candidateOrders.size()-1)
#endif
,_processDoneExternally(false)
{
    if(linearDetectors.size()!=_candidateOrders.size())
        throw RuntimeException("USIS::USIS: n� of detectors and number of channel matrix estimators (and candidate orders) are different.");

    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _linearDetectors[iChannelOrder] = linearDetectors[iChannelOrder]->Clone();
}


USIS::~USIS()
{
    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        delete _linearDetectors[iChannelOrder];

	delete _channelOrderEstimator;
}

vector<vector<tMatrix> > USIS::ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
	// channel estimation for the training sequence is needed in order to compute the channel order APP
	vector<vector<tMatrix> > estimatedMatrices = UnknownChannelOrderAlgorithm::ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

	// the estimated matrices are used to update the global channel order estimator
	tMatrix estimatedChannelOrderAPPs = _channelOrderEstimator->ComputeProbabilities(observations,estimatedMatrices,noiseVariances,trainingSequence);

	#ifdef CHANNELORDERSAPP_SAVING
		_channelOrderAPPs(_rCandidateOrders,tRange(_preamble.cols(),_preamble.cols()+trainingSequence.cols()-1)).inject(estimatedChannelOrderAPPs);
	#endif

    uint iChannelOrder;

    for(int i=_iFirstObservation;i<_iFirstObservation+trainingSequence.cols();i++)
    {
        for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        {
            // the observations from i to i+d are stacked
            tRange rSmoothingRange(i,i+_candidateOrders[iChannelOrder]-1);

            tVector stackedObservationsVector = Util::ToVector(observations(_rAllObservationRows,rSmoothingRange),columnwise);
            _linearDetectors[iChannelOrder]->StateStep(stackedObservationsVector);
        }

    }

    return estimatedMatrices;
}

void USIS::InitializeParticles()
{
    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter.Capacity();iParticle++)
    {
		// a clone of each of the channel matrix estimators...
		vector<ChannelMatrixEstimator *> thisParticleChannelMatrixEstimators(_candidateOrders.size());

		//...and linear detectors is constructed
		vector<LinearDetector *> thisParticleLinearDetectors(_candidateOrders.size());

		for(uint iCandidateOrder=0;iCandidateOrder<_candidateOrders.size();iCandidateOrder++)
		{
			thisParticleChannelMatrixEstimators[iCandidateOrder] = _channelEstimators[iCandidateOrder]->Clone();
			thisParticleLinearDetectors[iCandidateOrder] = _linearDetectors[iCandidateOrder]->Clone();
		}

		// ... and passed within a vector to each particle
		_particleFilter.AddParticle(new ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation(1.0/(double)_particleFilter.Capacity(),_N,_K,thisParticleChannelMatrixEstimators,thisParticleLinearDetectors,_channelOrderEstimator->Clone()));
    }
}

void USIS::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	int iParticle,iSmoothing,iRow,iSampledSymbol,iAlphabet,iSampled;
	uint iChannelOrder;
	int m,d,Nm,nLinearFiltersNeeded,iLinearFilterNeeded;
	vector<vector<tMatrix> > matricesToStack(_candidateOrders.size());
	tVector sampledVector(_N),sampledSmoothingVector(_N*_maxOrder);
	double proposal,s2q,sumProb,likelihoodsProd,sumLikelihoodsProd;
	vector<tMatrix> channelOrderEstimatorNeededSampledMatrices(_candidateOrders.size());

	// channel order APP saving
	int iBestParticle;

	// each matrix in "symbolProb" contains the probabilities connected to a channelOrder: symbolProb(i,j) is the p(i-th symbol=alphabet[j]). They are initialized with zeros
	vector<tMatrix> symbolProb(_candidateOrders.size(),LaGenMatDouble::zeros(_N*_maxOrder,_alphabet.Length()));

	// "overallSymbolProb" will combine the previous probabilities accordin to the APP of the channel order
	tMatrix overallSymbolProb(_N*_maxOrder,_alphabet.Length());

	// 2*_maxOrder-1 = m_{max} + d_{max}
	tMatrix forWeightUpdateNeededSymbols(_N,2*_maxOrder-1);

	// _maxOrder = d_{max} + 1
	tMatrix noiseCovariances[_maxOrder];

	tVector predictedNoiselessObservation(_L);

	int iObservationToBeProcessed = _startDetectionTime;
	while((iObservationToBeProcessed<_K) && !_processDoneExternally)
	{
		// observation matrix columns that are involved in the smoothing
		tRange rSmoothingRange(iObservationToBeProcessed,iObservationToBeProcessed+_maxOrder-1);

		// the stacked observations vector
		tVector stackedObservations = Util::ToVector(observations(_rAllObservationRows,rSmoothingRange),columnwise);

		// stacked noise covariance needs to be constructed
		tMatrix stackedNoiseCovariance = LaGenMatDouble::zeros(_L*_maxOrder,_L*_maxOrder);

		// the loop accomplishes 2 things:
		for(iSmoothing=0;iSmoothing<_maxOrder;iSmoothing++)
		{
			// i) construction of the stacked noise covariance
			for(iRow=0;iRow<_L;iRow++)
				stackedNoiseCovariance(iSmoothing*_L+iRow,iSmoothing*_L+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];

			// ii) obtaining the noise covariances for each time instant from the variances
			noiseCovariances[iSmoothing] = LaGenMatDouble::eye(_L);
			noiseCovariances[iSmoothing] *= noiseVariances[iObservationToBeProcessed+iSmoothing];
		}

		for(iParticle=0;iParticle<_particleFilter.Capacity();iParticle++)
		{
			ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *>(_particleFilter.GetParticle(iParticle));

			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				m = _candidateOrders[iChannelOrder];
				d = m-1;
				Nm = _N*m;
				matricesToStack[iChannelOrder] = vector<tMatrix>(_maxOrder,tMatrix(_L,Nm));

				tMatrix s2qAux(_L*(d+1),_L*(d+1));
				tVector s2qAuxFilter(_L*(d+1));

				// predicted channel matrices are sampled and stored in a vector in order to stack them

				// matricesToStack[iChannelOrder][0] = _ARcoefficient * <lastEstimatedChannelMatrix> + randn(_L,Nm)*_samplingVariance
				Util::Add(processedParticle->GetChannelMatrixEstimator(iChannelOrder)->LastEstimatedChannelMatrix(),StatUtil::RandnMatrix(_L,Nm,0.0,_samplingVariance),matricesToStack[iChannelOrder][0],_ARcoefficient,1.0);

				for(iSmoothing=1;iSmoothing<_maxOrder;iSmoothing++)
				{
					// matricesToStack[iChannelOrder][iSmoothing] = _ARcoefficient * matricesToStack[iChannelOrder][iSmoothing-1] + rand(_L,Nm)*_ARprocessVariance
					Util::Add(matricesToStack[iChannelOrder][iSmoothing-1],StatUtil::RandnMatrix(_L,Nm,0.0,_ARprocessVariance),matricesToStack[iChannelOrder][iSmoothing],_ARcoefficient,1.0);
				}

				// for sampling s_{t:t+d} we need to build
				nLinearFiltersNeeded = _maxOrder - m + 1; // linear filters

				// during the first iteration, we are going to use the real linear detector of this particle for this channel order
				LinearDetector *linearDetectorBeingProccessed = processedParticle->GetLinearDetector(iChannelOrder);

				for(iLinearFilterNeeded=0;iLinearFilterNeeded<nLinearFiltersNeeded;iLinearFilterNeeded++)
				{
					tRange rInvolvedObservations(iLinearFilterNeeded*_L,_L*(d+1+iLinearFilterNeeded)-1);

					// matrices are stacked to give
					tMatrix stackedChannelMatrix = HsToStackedH(matricesToStack[iChannelOrder],m,iLinearFilterNeeded,d+iLinearFilterNeeded);

					// the estimated stacked channel matrix is used to obtain soft estimations
					// of the transmitted symbols
					tVector softEstimations = linearDetectorBeingProccessed->Detect(stackedObservations(rInvolvedObservations),stackedChannelMatrix,stackedNoiseCovariance(rInvolvedObservations,rInvolvedObservations));

					tMatrix filter = linearDetectorBeingProccessed->ComputedFilter();

					// during the first iteration, we have used the real linear detector of this particle for this channel; during the remaining iterations we don't want the real linear detector to be modified
					if(iLinearFilterNeeded==0)
						linearDetectorBeingProccessed = linearDetectorBeingProccessed->Clone();

					// operations needed to computed the sampling variance

					//s2qAux = _alphabet.Variance() * stackedChannelMatrix * stackedChannelMatrix^H
            		Blas_Mat_Mat_Trans_Mult(stackedChannelMatrix,stackedChannelMatrix,s2qAux,_alphabet.Variance());

					// s2qAux = s2qAux + stackedNoiseCovariance
					Util::Add(s2qAux,stackedNoiseCovariance(rInvolvedObservations,rInvolvedObservations),s2qAux);

					// the real symbol we are sampling (it depends on "iLinearFilterNeeded")
					int iSampledSymbolPos = iLinearFilterNeeded*_N - 1;

					// sampling
					for(iSampledSymbol=0;iSampledSymbol<(_N*(d+1));iSampledSymbol++)
					{
						iSampledSymbolPos++;

						// s2qAuxFilter = s2qAux * filter.col(iSampledSymbol)
						Blas_Mat_Vec_Mult(s2qAux,filter.col(iSampledSymbol),s2qAuxFilter);

                		s2q = _alphabet.Variance()*(1.0 - 2.0*Blas_Dot_Prod(filter.col(iSampledSymbol),stackedChannelMatrix.col(_N*(m-1)+iSampledSymbol))) + Blas_Dot_Prod(filter.col(iSampledSymbol),s2qAuxFilter);

						sumProb = 0.0;
						// the probability for each posible symbol alphabet is computed
						for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
						{
							symbolProb[iChannelOrder](iSampledSymbolPos,iAlphabet) = StatUtil::NormalPdf(softEstimations(iSampledSymbol),_alphabet[iAlphabet],s2q);

							// the computed pdf is accumulated for normalizing purposes
							sumProb += symbolProb[iChannelOrder](iSampledSymbolPos,iAlphabet);
						}

						try {
							for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
								symbolProb[iChannelOrder](iSampledSymbolPos,iAlphabet) /= sumProb;
						}catch(exception e){
							cout << "The sum of the probabilities is null." << endl;
							for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
								symbolProb[iChannelOrder](iSampledSymbolPos,iAlphabet) = 0.5;
						}
					}
				} // for(iLinearFilterNeeded=0;iLinearFilterNeeded<nLinearFiltersNeeded;iLinearFilterNeeded++)

				// the clone is dismissed
				delete linearDetectorBeingProccessed;

			} //for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)

			//the probabilities of the different channel orders are weighted according to the a posteriori probability of the channel order in the previous time instant
			overallSymbolProb = symbolProb[0];
			overallSymbolProb *= processedParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(0);
			for(iChannelOrder=1;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				Util::Add(overallSymbolProb,symbolProb[iChannelOrder],overallSymbolProb,1.0,processedParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(iChannelOrder));
			}

			proposal = 1.0;
			// the symbols are sampled from the above combined probabilities
			for(iSampledSymbol=0;iSampledSymbol<_N*_maxOrder;iSampledSymbol++)
			{
				iSampled = StatUtil::Discrete_rnd(overallSymbolProb.row(iSampledSymbol));

				sampledSmoothingVector(iSampledSymbol) = _alphabet[iSampled];

				proposal *= overallSymbolProb(iSampledSymbol,iSampled);
			}

			// sampled symbol vector is stored for the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledSmoothingVector(_allSymbolsRows));

			sumLikelihoodsProd = 0.0;

			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				m = _candidateOrders[iChannelOrder];
				d = m-1;

				// range for the already detected symbol vectors involved in the current detection
				tRange rAlreadyDetectedSymbolVectors(iObservationToBeProcessed-m+1,iObservationToBeProcessed-1);

				tRange r0mMinus2(0,m-2),rSampledSymbolVectors(m-1,m+_maxOrder-2);
				tRange rFirstmSymbolVectors(0,m-1);

				// all the symbol vectors involved in the smoothing are kept in "forWeightUpdateNeededSymbols":
				// i) the already known:
				forWeightUpdateNeededSymbols(_allSymbolsRows,r0mMinus2).inject(processedParticle->GetSymbolVectors(rAlreadyDetectedSymbolVectors));

				// ii) the just sampled
				forWeightUpdateNeededSymbols(_allSymbolsRows,rSampledSymbolVectors).inject(Util::ToMatrix(sampledSmoothingVector,columnwise,_N));

				likelihoodsProd = processedParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(iChannelOrder);

				for(iSmoothing=0;iSmoothing<_maxOrder;iSmoothing++)
				{
					tRange rSymbolVectors(iSmoothing,iSmoothing+m-1);
					tVector stackedSymbolVector = Util::ToVector(forWeightUpdateNeededSymbols(_allSymbolsRows,rSymbolVectors),columnwise);

					// predictedNoiselessObservation = matricesToStack[iChannelOrder][iSmoothing] * stackedSymbolVector
					Blas_Mat_Vec_Mult(matricesToStack[iChannelOrder][iSmoothing],stackedSymbolVector,predictedNoiselessObservation);

					likelihoodsProd *= StatUtil::NormalPdf(observations.col(iObservationToBeProcessed+iSmoothing),predictedNoiselessObservation,noiseCovariances[iSmoothing]);

				} // for(iSmoothing=0;iSmoothing<_maxOrder;iSmoothing++)

				// the estimation of the channel matrix is updated
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed, (processedParticle->GetChannelMatrixEstimator(iChannelOrder))->NextMatrix(observations.col(iObservationToBeProcessed),forWeightUpdateNeededSymbols(_allSymbolsRows,rFirstmSymbolVectors),noiseVariances[iObservationToBeProcessed]));

                // the computed likelihood is accumulated
                sumLikelihoodsProd += likelihoodsProd;

				// a vector of channel matrices is built to update the channel order estimator
                channelOrderEstimatorNeededSampledMatrices[iChannelOrder] = matricesToStack[iChannelOrder][0];
			} // for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)

			// the channel order estimator is updated
			processedParticle->GetChannelOrderEstimator()->Update(observations.col(iObservationToBeProcessed),channelOrderEstimatorNeededSampledMatrices,sampledSmoothingVector(_allSymbolsRows),noiseVariances[iObservationToBeProcessed]);

			// the weight is updated
			processedParticle->SetWeight((sumLikelihoodsProd/proposal)*processedParticle->GetWeight());

		} // for(iParticle=0;iParticle<_particleFilter.Capacity();iParticle++)

		_particleFilter.NormalizeWeights();

		// we find out which is the "best" particle at this time instant
		Util::Max(_particleFilter.GetWeightsVector(),iBestParticle);

		#ifdef CHANNELORDERSAPP_SAVING
			// its a posteriori channel order probabilities are stored
			ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *bestParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *>(_particleFilter.GetParticle(iBestParticle));

			for(uint i=0;i<_candidateOrders.size();i++)
				_channelOrderAPPs(i,iObservationToBeProcessed) = bestParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(i);
		#endif

		BeforeResamplingProcess(iObservationToBeProcessed,observations,noiseVariances);

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
            _resamplingAlgorithm->ResampleWhenNecessary(&_particleFilter);

    	iObservationToBeProcessed++;
	} // while((iObservationToBeProcessed<_K) && !_processDoneExternally)
}


int USIS::BestChannelOrderIndex(int iBestParticle)
{
	ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *bestParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *>(_particleFilter.GetParticle(iBestParticle));

	int iMaxChannelOrderAPP = 0;
	double maxChannelOrderAPP = bestParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(iMaxChannelOrderAPP);

	for(uint i=1;i<_candidateOrders.size();i++)
		if(bestParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(i) > maxChannelOrderAPP)
		{
			maxChannelOrderAPP = bestParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(i);
			iMaxChannelOrderAPP = i;
		}

	return iMaxChannelOrderAPP;
}
