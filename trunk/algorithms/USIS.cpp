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

// #define DEBUG12

USIS::USIS(string name, Alphabet alphabet, int L, int N, int frameLength, vector< ChannelMatrixEstimator * > channelEstimators,vector<LinearDetector *> linearDetectors, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,ChannelOrderEstimator * channelOrderEstimator,double ARcoefficient,double samplingVariance,double ARprocessVariance): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, frameLength, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_linearDetectors(linearDetectors.size()),_channelOrderEstimator(channelOrderEstimator->Clone()),_particleFilter(nParticles),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance),_rAllObservationRows(0,_nOutputs-1)
,_processDoneExternally(false)
{
    if(linearDetectors.size()!=_candidateOrders.size())
        throw RuntimeException("USIS::USIS: number of detectors and number of channel matrix estimators (and candidate orders) are different.");

    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _linearDetectors[iChannelOrder] = linearDetectors[iChannelOrder]->Clone();
}


USIS::~USIS()
{
    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        delete _linearDetectors[iChannelOrder];

	delete _channelOrderEstimator;
}

// vector<vector<tMatrix> > USIS::EstimateChannelFromTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
// {
// 	// channel estimation for the training sequence is needed in order to compute the channel order APP
// 	vector<vector<tMatrix> > estimatedMatrices = UnknownChannelOrderAlgorithm::EstimateChannelFromTrainingSequence(observations,noiseVariances,trainingSequence);
//
// 	// the estimated matrices are used to update the global channel order estimator and compute the channel order APP
//     // during the training sequence
// 	tMatrix estimatedChannelOrderAPPs = _channelOrderEstimator->ComputeProbabilities(observations,estimatedMatrices,noiseVariances,Util::append(_preamble,trainingSequence),_preamble.cols());
//
//     // the APP of the candidate channel orders are set accordingly
// 	_channelOrderAPPs(tRange(),tRange(_preamble.cols(),_preamble.cols()+trainingSequence.cols()-1)).inject(estimatedChannelOrderAPPs);
//
//     return estimatedMatrices;
// }

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

            if(_randomParticlesInitilization)
                // the first matrix of the channel matrix estimator is initialized randomly
                thisParticleChannelMatrixEstimators[iCandidateOrder]->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::RandMatrix(_channelMeanVectors[iCandidateOrder],_channelCovariances[iCandidateOrder],StatUtil::_particlesInitializerRandomGenerator),rowwise,_nOutputs));

			thisParticleLinearDetectors[iCandidateOrder] = _linearDetectors[iCandidateOrder]->Clone();
		}

		// ... and passed within a vector to each particle
		_particleFilter.AddParticle(new ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation(1.0/(double)_particleFilter.Capacity(),_nInputs,_K,thisParticleChannelMatrixEstimators,thisParticleLinearDetectors,_channelOrderEstimator->Clone()));
    }
}

void USIS::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	int iParticle,iSmoothing,iRow,iSampledSymbol,iAlphabet,iSampled;
	uint iChannelOrder;
	int m,d,Nm,nLinearFiltersNeeded,iLinearFilterNeeded;
	vector<vector<tMatrix> > matricesToStack(_candidateOrders.size());
	tVector sampledVector(_nInputs),sampledSmoothingVector(_nInputs*_maxOrder);
	double proposal,s2q,sumProb,likelihoodsProd,sumLikelihoodsProd;
	vector<tMatrix> channelOrderEstimatorNeededSampledMatrices(_candidateOrders.size());

	// each matrix in "symbolProb" contains the probabilities connected to a channelOrder: symbolProb(i,j) is the p(i-th symbol=alphabet[j]). They are initialized with zeros
	vector<tMatrix> symbolProb(_candidateOrders.size(),LaGenMatDouble::zeros(_nInputs*_maxOrder,_alphabet.length()));

	// "overallSymbolProb" will combine the previous probabilities accordin to the APP of the channel order
	tMatrix overallSymbolProb(_nInputs*_maxOrder,_alphabet.length());

	// 2*_maxOrder-1 = m_{max} + d_{max}
	tMatrix forWeightUpdateNeededSymbols(_nInputs,2*_maxOrder-1);

	// _maxOrder = d_{max} + 1
	tMatrix noiseCovariances[_maxOrder];

	tVector predictedNoiselessObservation(_nOutputs);

	int iObservationToBeProcessed = _startDetectionTime;
	while((iObservationToBeProcessed<_K) && !_processDoneExternally)
	{
#ifdef DEBUG
		cout << "iObservationToBeProcessed = " << iObservationToBeProcessed << endl;
		cout << "_K = " << _K << endl;

		if(iObservationToBeProcessed>20)
			exit(0);
#endif

		// observation matrix columns that are involved in the smoothing
		tRange rSmoothingRange(iObservationToBeProcessed,iObservationToBeProcessed+_maxOrder-1);

		// the stacked observations vector
		tVector stackedObservations = Util::toVector(observations(_rAllObservationRows,rSmoothingRange),columnwise);

		// stacked noise covariance needs to be constructed
		tMatrix stackedNoiseCovariance = LaGenMatDouble::zeros(_nOutputs*_maxOrder,_nOutputs*_maxOrder);

		// the loop accomplishes 2 things:
		for(iSmoothing=0;iSmoothing<_maxOrder;iSmoothing++)
		{
			// i) construction of the stacked noise covariance
			for(iRow=0;iRow<_nOutputs;iRow++)
				stackedNoiseCovariance(iSmoothing*_nOutputs+iRow,iSmoothing*_nOutputs+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];

			// ii) obtaining the noise covariances for each time instant from the variances
			noiseCovariances[iSmoothing] = LaGenMatDouble::eye(_nOutputs);
			noiseCovariances[iSmoothing] *= noiseVariances[iObservationToBeProcessed+iSmoothing];
		}

		for(iParticle=0;iParticle<_particleFilter.Capacity();iParticle++)
		{
#ifdef DEBUG12
            cout << "iParticle = " << iParticle << endl;
#endif
			ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *>(_particleFilter.GetParticle(iParticle));

			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				m = _candidateOrders[iChannelOrder];
				d = m-1;
				Nm = _nInputs*m;
				matricesToStack[iChannelOrder] = vector<tMatrix>(_maxOrder,tMatrix(_nOutputs,Nm));

				tMatrix s2qAux(_nOutputs*(d+1),_nOutputs*(d+1));
				tVector s2qAuxFilter(_nOutputs*(d+1));

				// predicted channel matrices are sampled and stored in a vector in order to stack them

				// matricesToStack[iChannelOrder][0] = _ARcoefficient * <lastEstimatedChannelMatrix> + randn(_nOutputs,Nm)*_samplingVariance
				Util::add(processedParticle->GetChannelMatrixEstimator(iChannelOrder)->lastEstimatedChannelMatrix(),StatUtil::RandnMatrix(_nOutputs,Nm,0.0,_samplingVariance),matricesToStack[iChannelOrder][0],_ARcoefficient,1.0);

				for(iSmoothing=1;iSmoothing<_maxOrder;iSmoothing++)
				{
					// matricesToStack[iChannelOrder][iSmoothing] = _ARcoefficient * matricesToStack[iChannelOrder][iSmoothing-1] + rand(_nOutputs,Nm)*_ARprocessVariance
					Util::add(matricesToStack[iChannelOrder][iSmoothing-1],StatUtil::RandnMatrix(_nOutputs,Nm,0.0,_ARprocessVariance),matricesToStack[iChannelOrder][iSmoothing],_ARcoefficient,1.0);
				}

				// for sampling s_{t:t+d} we need to build
				nLinearFiltersNeeded = _maxOrder - m + 1; // linear filters

				// during the first iteration, we are going to use the real linear detector of this particle for this channel order
				LinearDetector *linearDetectorBeingProccessed = processedParticle->GetLinearDetector(iChannelOrder);

				for(iLinearFilterNeeded=0;iLinearFilterNeeded<nLinearFiltersNeeded;iLinearFilterNeeded++)
				{
					tRange rInvolvedObservations(iLinearFilterNeeded*_nOutputs,_nOutputs*(d+1+iLinearFilterNeeded)-1);

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

					//s2qAux = _alphabet.variance() * stackedChannelMatrix * stackedChannelMatrix^H
            		Blas_Mat_Mat_Trans_Mult(stackedChannelMatrix,stackedChannelMatrix,s2qAux,_alphabet.variance());

					// s2qAux = s2qAux + stackedNoiseCovariance
					Util::add(s2qAux,stackedNoiseCovariance(rInvolvedObservations,rInvolvedObservations),s2qAux);

					// the real symbol we are sampling (it depends on "iLinearFilterNeeded")
					int iSampledSymbolPos = iLinearFilterNeeded*_nInputs - 1;

					// sampling
					for(iSampledSymbol=0;iSampledSymbol<(_nInputs*(d+1));iSampledSymbol++)
					{
						iSampledSymbolPos++;

						// s2qAuxFilter = s2qAux * filter.col(iSampledSymbol)
						Blas_Mat_Vec_Mult(s2qAux,filter.col(iSampledSymbol),s2qAuxFilter);

                		s2q = _alphabet.variance()*(1.0 - 2.0*Blas_Dot_Prod(filter.col(iSampledSymbol),stackedChannelMatrix.col(_nInputs*(m-1)+iSampledSymbol))) + Blas_Dot_Prod(filter.col(iSampledSymbol),s2qAuxFilter);

						sumProb = 0.0;
						// the probability for each posible symbol alphabet is computed
						for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
						{
							symbolProb[iChannelOrder](iSampledSymbolPos,iAlphabet) = StatUtil::NormalPdf(softEstimations(iSampledSymbol),_alphabet[iAlphabet],s2q);

							// the computed pdf is accumulated for normalizing purposes
							sumProb += symbolProb[iChannelOrder](iSampledSymbolPos,iAlphabet);
						}

						try {
							for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
								symbolProb[iChannelOrder](iSampledSymbolPos,iAlphabet) /= sumProb;
						}catch(exception e){
							cout << "The sum of the probabilities is null." << endl;
							for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
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
				Util::add(overallSymbolProb,symbolProb[iChannelOrder],overallSymbolProb,1.0,processedParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(iChannelOrder));
			}

			proposal = 1.0;
			// the symbols are sampled from the above combined probabilities
			for(iSampledSymbol=0;iSampledSymbol<_nInputs*_maxOrder;iSampledSymbol++)
			{
				iSampled = StatUtil::discrete_rnd(overallSymbolProb.row(iSampledSymbol));

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
				forWeightUpdateNeededSymbols(_allSymbolsRows,rSampledSymbolVectors).inject(Util::toMatrix(sampledSmoothingVector,columnwise,_nInputs));

				likelihoodsProd = processedParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(iChannelOrder);

				for(iSmoothing=0;iSmoothing<_maxOrder;iSmoothing++)
				{
					tRange rSymbolVectors(iSmoothing,iSmoothing+m-1);
					tVector stackedSymbolVector = Util::toVector(forWeightUpdateNeededSymbols(_allSymbolsRows,rSymbolVectors),columnwise);

					// predictedNoiselessObservation = matricesToStack[iChannelOrder][iSmoothing] * stackedSymbolVector
					Blas_Mat_Vec_Mult(matricesToStack[iChannelOrder][iSmoothing],stackedSymbolVector,predictedNoiselessObservation);

					likelihoodsProd *= StatUtil::NormalPdf(observations.col(iObservationToBeProcessed+iSmoothing),predictedNoiselessObservation,noiseCovariances[iSmoothing]);

				} // for(iSmoothing=0;iSmoothing<_maxOrder;iSmoothing++)

				// the estimation of the channel matrix is updated
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed,
				(processedParticle->GetChannelMatrixEstimator(iChannelOrder))->nextMatrix(observations.col(iObservationToBeProcessed),
					forWeightUpdateNeededSymbols(_allSymbolsRows,rFirstmSymbolVectors),noiseVariances[iObservationToBeProcessed]));

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
		ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *bestParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *>(_particleFilter.GetBestParticle());

        // its a posteriori channel order probabilities are stored
		for(uint i=0;i<_candidateOrders.size();i++)
			_channelOrderAPPs(i,iObservationToBeProcessed) = bestParticle->GetChannelOrderEstimator()->GetChannelOrderAPP(i);

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

void USIS::BeforeInitializingParticles(const tMatrix &observations,vector<double> &noiseVariances,const tMatrix &trainingSequence)
{
//     for(int i=_iFirstObservation;i<_iFirstObservation+trainingSequence.cols();i++)
//     {
//         for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
//         {
//             // the observations from i to i+d are stacked
//             tRange rSmoothingRange(i,i+_candidateOrders[iChannelOrder]-1);
//
//             tVector stackedObservationsVector = Util::toVector(observations(_rAllObservationRows,rSmoothingRange),columnwise);
//             _linearDetectors[iChannelOrder]->StateStep(stackedObservationsVector);
//         }
//
//     }

    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _linearDetectors[iChannelOrder]->StateStepsFromObservationsSequence(observations,_candidateOrders[iChannelOrder]-1,_preamble.cols(),_preamble.cols()+trainingSequence.cols());
//     {
//         // the observations from i to i+d are stacked
//         tRange rSmoothingRange(i,i+_candidateOrders[iChannelOrder]-1);
//
//         tVector stackedObservationsVector = Util::toVector(observations(_rAllObservationRows,rSmoothingRange),columnwise);
//         _linearDetectors[iChannelOrder]->StateStep(stackedObservationsVector);
//     }

    // the APP of the candidate channel orders are set accordingly
    _channelOrderAPPs(tRange(),tRange(_preamble.cols(),_preamble.cols()+trainingSequence.cols()-1)) = 1.0/double(_candidateOrders.size());
}

void USIS::UpdateParticleChannelOrderEstimators(Particle *particle,const tMatrix &observations,const std::vector<std::vector<tMatrix> > &channelMatrices,vector<double> &noiseVariances,const tMatrix &sequenceToProcess)
{
    ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation* convertedParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *> (particle);

    convertedParticle->GetChannelOrderEstimator()->ComputeProbabilities(observations,channelMatrices,noiseVariances,sequenceToProcess,_preamble.cols());
}
