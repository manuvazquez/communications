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
#include "LinearFilterBasedUnknownChannelOrderSMCAlgorithm.h"

// #define DEBUG4

LinearFilterBasedUnknownChannelOrderSMCAlgorithm::LinearFilterBasedUnknownChannelOrderSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators,vector<LinearDetector *> linearDetectors, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance/*,const MIMOChannel &canal,const tMatrix &simbolos*/): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_allObservationRows(0,_L-1),_linearDetectors(linearDetectors.size()),_particleFilter(nParticles),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance)
// ,_canal(canal),_simbolos(simbolos)
{
    if(linearDetectors.size()!=_candidateOrders.size())
        throw RuntimeException("LinearFilterBasedUnknownChannelOrderSMCAlgorithm::LinearFilterBasedUnknownChannelOrderSMCAlgorithm: nº of detectors and number of channel matrix estimators (and candidate orders) are different.");

    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _linearDetectors[iChannelOrder] = linearDetectors[iChannelOrder]->Clone();
}


LinearFilterBasedUnknownChannelOrderSMCAlgorithm::~LinearFilterBasedUnknownChannelOrderSMCAlgorithm()
{
    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        delete _linearDetectors[iChannelOrder];
}

vector<vector<tMatrix> > LinearFilterBasedUnknownChannelOrderSMCAlgorithm::ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
    int iChannelOrder;

    for(int i=_iFirstObservation;i<_iFirstObservation+trainingSequence.cols();i++)
    {
        for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        {
            // the observations from i to i+d are stacked
            tRange rSmoothingRange(i,i+_candidateOrders[iChannelOrder]-1);
            tVector stackedObservationsVector = Util::ToVector(observations(_allObservationRows,rSmoothingRange),columnwise);
            _linearDetectors[iChannelOrder]->StateStep(stackedObservationsVector);
        }
    }

    return UnknownChannelOrderAlgorithm::ProcessTrainingSequence(observations,noiseVariances,trainingSequence);
}

void LinearFilterBasedUnknownChannelOrderSMCAlgorithm::InitializeParticles()
{
    // memory is reserved
    for(int iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
    {
		// a clone of each of the channel matrix estimators...
		vector<ChannelMatrixEstimator *> thisParticleChannelMatrixEstimators(_candidateOrders.size());

		//...and linear detectors is constructed
		vector<LinearDetector *> thisParticleLinearDetectors(_candidateOrders.size());

		for(int iCandidateOrder=0;iCandidateOrder<_candidateOrders.size();iCandidateOrder++)
		{
			thisParticleChannelMatrixEstimators[iCandidateOrder] = _channelEstimators[iCandidateOrder]->Clone();
			thisParticleLinearDetectors[iCandidateOrder] = _linearDetectors[iCandidateOrder]->Clone();
		}

		// ... and passed within a vector to each particle
		_particleFilter.SetParticle(new ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP(1.0/(double)_particleFilter.Nparticles(),_N,_K,thisParticleChannelMatrixEstimators,_candidateOrders.size(),thisParticleLinearDetectors),iParticle);
    }
}

void LinearFilterBasedUnknownChannelOrderSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	int iParticle,iSmoothing,iRow,iSampledSymbol,iAlphabet,iSampled,iChannelOrder;
	int m,d,Nm,nLinearFiltersNeeded,iLinearFilterNeeded;
	vector<vector<tMatrix> > matricesToStack(_candidateOrders.size());
	tVector sampledVector(_N),sampledSmoothingVector(_N*_maxOrder);
	double proposal,s2q,sumProb,likelihoodsProd,sumLikelihoodsProd,channelOrderAPPsNormConstant;
	double *newChannelOrderAPPs = new double[_candidateOrders.size()];

	// each matrix in "symbolProb" contains the probabilities connected to a channelOrder: symbolProb(i,j) is the p(i-th symbol=alphabet[j]). They are initialized with zeros
	tMatrix symbolProbAux = LaGenMatDouble::zeros(_N*_maxOrder,_alphabet.Length());
	vector<tMatrix> symbolProb(_candidateOrders.size(),symbolProbAux);

	// "overallSymbolProb" will combine the previous probabilities accordin to the APP of the channel order
	tMatrix overallSymbolProb(_N*_maxOrder,_alphabet.Length());

	// 2*_maxOrder-1 = m_{max} + d_{max}
	tMatrix forWeightUpdateNeededSymbols(_N,2*_maxOrder-1);

	tMatrix noiseCovariances[_maxOrder];
	tVector predictedNoiselessObservation(_L);

	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
		#ifdef DEBUG2
			cout << "iObservationToBeProcessed: " << iObservationToBeProcessed << endl;
			cout << "Simbolos anteriores" << endl << _simbolos(_allSymbolsRows, tRange(iObservationToBeProcessed-_maxOrder+1,iObservationToBeProcessed-1)) << endl;
			cout << "Simbolos transmitidos" << endl << _simbolos(_allSymbolsRows, tRange(iObservationToBeProcessed,iObservationToBeProcessed+_maxOrder-1)) << endl;
		#endif

		// observation matrix columns that are involved in the smoothing
		tRange rSmoothingRange(iObservationToBeProcessed,iObservationToBeProcessed+_maxOrder-1);

		// the stacked observations vector
		tVector stackedObservations = Util::ToVector(observations(_allObservationRows,rSmoothingRange),columnwise);

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

		for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
		{
			ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP *>(_particleFilter.GetParticle(iParticle));

			#ifdef DEBUG
				for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
				{
					cout << "Orden " << _candidateOrders[iChannelOrder] << " " << processedParticle->GetChannelOrderAPP(iChannelOrder) << " ";
				}
				cout << endl;
			#endif

			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				m = _candidateOrders[iChannelOrder];
				d = m-1;
				Nm = _N*m;
				matricesToStack[iChannelOrder] = vector<tMatrix>(_maxOrder,tMatrix(_L,Nm));

				#ifdef DEBUG4
					cout << "m = " << m << endl;
				#endif

				tMatrix s2qAux(_L*(d+1),_L*(d+1));
				tVector s2qAuxFilter(_L*(d+1));

				// predicted channel matrices are sampled and stored in a vector in order to stack them
				// matricesToStack[iChannelOrder][0] = _ARcoefficient * <lastEstimatedChannelMatrix> + randn(_L,Nm)*_samplingVariance
				Util::Add((processedParticle->GetChannelMatrixEstimator(iChannelOrder))->LastEstimatedChannelMatrix(),StatUtil::RandnMatrix(_L,Nm,0.0,_samplingVariance),matricesToStack[iChannelOrder][0],_ARcoefficient,1.0);

				for(iSmoothing=1;iSmoothing<_maxOrder;iSmoothing++)
				{
					// matricesToStack[iChannelOrder][iSmoothing] = _ARcoefficient * matricesToStack[iChannelOrder][iSmoothing-1] + rand(_L,Nm)*_ARprocessVariance
					Util::Add(matricesToStack[iChannelOrder][iSmoothing-1],StatUtil::RandnMatrix(_L,Nm,0.0,_ARprocessVariance),matricesToStack[iChannelOrder][iSmoothing],_ARcoefficient,1.0);
				}

				// for sampling s_{t:t+d} we need to build
				nLinearFiltersNeeded = _maxOrder - m + 1; // linear filters

// 				nLinearFiltersNeeded = 1;

				// during the first iteration, we are going to use the real linear detector of this particle for this channel order
				LinearDetector *linearDetectorBeingProccessed = processedParticle->GetLinearDetector(iChannelOrder);

				#ifdef DEBUG3
					cout << "Simbolos" << endl << _simbolos(_allSymbolsRows, tRange(iObservationToBeProcessed,iObservationToBeProcessed+_maxOrder-1)) << endl;
				#endif

				for(iLinearFilterNeeded=0;iLinearFilterNeeded<nLinearFiltersNeeded;iLinearFilterNeeded++)
				{
					#ifdef DEBUG4
						cout << "iLinearFilterNeeded = " << iLinearFilterNeeded << endl;
					#endif

					tRange rInvolvedObservations(iLinearFilterNeeded*_L,_L*(d+1+iLinearFilterNeeded)-1);
					// matrices are stacked to give
					tMatrix stackedChannelMatrix = HsToStackedH(matricesToStack[iChannelOrder],m,iLinearFilterNeeded,d+iLinearFilterNeeded);

					#ifdef DEBUG4
						cout << "stackedChannelMatrix" << endl << stackedChannelMatrix << endl;
						cout << "stackedNoiseCovariance(rInvolvedObservations,rInvolvedObservations)" << endl << stackedNoiseCovariance(rInvolvedObservations,rInvolvedObservations) << endl;
					#endif

//             		tMatrix stackedChannelMatrixMinus = stackedChannelMatrix(tRange(0,_L*(d+1)-1),tRange((m-1)*_N,stackedChannelMatrix.cols()-1));

					// the estimated stacked channel matrix is used to obtain soft estimations
					// of the transmitted symbols
					tVector softEstimations = linearDetectorBeingProccessed->Detect(stackedObservations(rInvolvedObservations),stackedChannelMatrix,stackedNoiseCovariance(rInvolvedObservations,rInvolvedObservations));

					#ifdef DEBUG3
						cout << "iLinearFilterNeeded = " << iLinearFilterNeeded << endl << softEstimations << endl;
						getchar();
					#endif

					tMatrix filter = linearDetectorBeingProccessed->ComputedFilter();

					#ifdef DEBUG4
						cout << "filter" << endl << filter << endl;
					#endif

					// during the first iteration, we have used the real linear detector of this particle for this channel; during the remaining iterations we don't want the real linear detector to be modified
					if(iLinearFilterNeeded==0)
						linearDetectorBeingProccessed = linearDetectorBeingProccessed->Clone();

					// operations needed to computed the sampling variance

					//s2qAux = _alphabet.Variance() * stackedChannelMatrix * stackedChannelMatrix^H
            		Blas_Mat_Mat_Trans_Mult(stackedChannelMatrix,stackedChannelMatrix,s2qAux,_alphabet.Variance());
//             		Blas_Mat_Mat_Trans_Mult(stackedChannelMatrixMinus,stackedChannelMatrixMinus,s2qAux,_alphabet.Variance());

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

//                 		s2q = _alphabet.Variance()*(1.0 - 2.0*Blas_Dot_Prod(filter.col(iSampledSymbol),stackedChannelMatrixMinus.col(iSampledSymbol))) + Blas_Dot_Prod(filter.col(iSampledSymbol),s2qAuxFilter);

// 						s2q /= 10;

						#ifdef DEBUG4
							cout << "s2qAuxFilter" << endl << s2qAuxFilter << endl;
                			cout << "s2q = " << s2q << endl;
                			getchar();
                		#endif

						double sumProb = 0.0;
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

				#ifdef DEBUG4
					cout << "probabilidades par orden " << _candidateOrders[iChannelOrder] << endl << symbolProb[iChannelOrder] << endl;
				#endif

				// the clone is dismissed
				delete linearDetectorBeingProccessed;

				#ifdef DEBUG4
					getchar();
				#endif

			} //for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)

			//the probabilities of the different channel orders are weighted according to the a posteriori probability of the channel order in the previous time instant
			overallSymbolProb = symbolProb[0];
			overallSymbolProb *= processedParticle->GetChannelOrderAPP(0);
			for(iChannelOrder=1;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				Util::Add(overallSymbolProb,symbolProb[iChannelOrder],overallSymbolProb,1.0,processedParticle->GetChannelOrderAPP(iChannelOrder));
			}

			proposal = 1.0;
			// the symbols are sampled from the above combined probabilities
// 			for(iSampledSymbol=0;iSampledSymbol<_N;iSampledSymbol++)
			for(iSampledSymbol=0;iSampledSymbol<_N*_maxOrder;iSampledSymbol++)
			{
				int iSampled = StatUtil::Discrete_rnd(overallSymbolProb.row(iSampledSymbol));

// 				cout << "iSampled fue " << iSampled << endl;
// 				getchar();
				sampledSmoothingVector(iSampledSymbol) = _alphabet[iSampled];

				proposal *= overallSymbolProb(iSampledSymbol,iSampled);
			}

			#ifdef DEBUG
				for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
				{
					cout << " Orden " << _candidateOrders[iChannelOrder] << endl << symbolProb[iChannelOrder] << endl;
				}
				cout << "Se va a nuestrear de" << endl << overallSymbolProb << endl;
			#endif

			// sampled symbol vector is stored for the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledSmoothingVector(_allSymbolsRows));

			#ifdef DEBUG
				cout << "Muestreó" << endl << processedParticle->GetSymbolVector(iObservationToBeProcessed) << endl;
			#endif

			sumLikelihoodsProd = 0.0;
			channelOrderAPPsNormConstant = 0.0;

			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				m = _candidateOrders[iChannelOrder];
				d = m-1;

				// range for the already detected symbol vectors involved in the current detection
				tRange rAlreadyDetectedSymbolVectors(iObservationToBeProcessed-m+1,iObservationToBeProcessed-1);

				tRange r0mMinus2(0,m-2),rSampledSymbolVectors(m-1,m+_maxOrder-2);
				tRange rFirstmSymbolVectors(0,m-1);

				// all the symbol vectors involved in the smoothing are kept in "forWeightUpdateNeededSymbols"
				// i) the already known:
				forWeightUpdateNeededSymbols(_allSymbolsRows,r0mMinus2).inject(processedParticle->GetSymbolVectors(rAlreadyDetectedSymbolVectors));

				// ii) the just sampled
				forWeightUpdateNeededSymbols(_allSymbolsRows,rSampledSymbolVectors).inject(Util::ToMatrix(sampledSmoothingVector,columnwise,_N));

				likelihoodsProd = processedParticle->GetChannelOrderAPP(iChannelOrder);

				#ifdef DEBUG3
					cout << "La anterior probabilidad de m para orden " << _candidateOrders[iChannelOrder] << ": " << likelihoodsProd << endl;
				#endif

// 				for(iSmoothing=0;iSmoothing<1;iSmoothing++)
				for(iSmoothing=0;iSmoothing<_maxOrder;iSmoothing++)
				{
					tRange rSymbolVectors(iSmoothing,iSmoothing+m-1);
					tVector stackedSymbolVector = Util::ToVector(forWeightUpdateNeededSymbols(_allSymbolsRows,rSymbolVectors),columnwise);

					#ifdef DEBUG2
						cout << "iSmoothing = " << iSmoothing << endl;
						cout << "Simbolos para verosimilitud: " << endl << forWeightUpdateNeededSymbols(_allSymbolsRows,rSymbolVectors);
						cout << "Matriz de canal: " << endl << matricesToStack[iChannelOrder][iSmoothing];
					#endif

					// predictedNoiselessObservation = matricesToStack[iChannelOrder][iSmoothing] * stackedSymbolVector
					Blas_Mat_Vec_Mult(matricesToStack[iChannelOrder][iSmoothing],stackedSymbolVector,predictedNoiselessObservation);

					likelihoodsProd *= StatUtil::NormalPdf(observations.col(iObservationToBeProcessed+iSmoothing),predictedNoiselessObservation,noiseCovariances[iSmoothing]);

					// unnormalized APP channel order
					if(iSmoothing==0)
					{
// 						processedParticle->SetChannelOrderAPP(likelihoodsProd,iChannelOrder);
						newChannelOrderAPPs[iChannelOrder] = likelihoodsProd;

						// it's accumulated for normalization purposes
						channelOrderAPPsNormConstant += likelihoodsProd;

						#ifdef DEBUG3
							cout << "Se actualiza la probabilidad de m por " << StatUtil::NormalPdf(observations.col(iObservationToBeProcessed+iSmoothing),predictedNoiselessObservation,noiseCovariances[iSmoothing]) << endl;
						#endif
					}
				}

				// the estimation of the channel matrix is updated
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed, (processedParticle->GetChannelMatrixEstimator(iChannelOrder))->NextMatrix(observations.col(iObservationToBeProcessed),forWeightUpdateNeededSymbols(_allSymbolsRows,rFirstmSymbolVectors),noiseVariances[iObservationToBeProcessed]));

                // the computed likelihood is accumulated
                sumLikelihoodsProd += likelihoodsProd;
			}

			// the channel order APPs are normalized for the next iteration
			if(channelOrderAPPsNormConstant!=0)
			{
				for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
					processedParticle->SetChannelOrderAPP(newChannelOrderAPPs[iChannelOrder]/channelOrderAPPsNormConstant,iChannelOrder);
			}

// 			// the channel order APPs are normalized for the next iteration
// 			if(channelOrderAPPsNormConstant==0)
// 			{
// 				for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
// 					processedParticle->SetChannelOrderAPP(1.0/(double)_candidateOrders.size(),iChannelOrder);
// 			}
// 			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
// 			{
// 				processedParticle->SetChannelOrderAPP(processedParticle->GetChannelOrderAPP(iChannelOrder)/channelOrderAPPsNormConstant,iChannelOrder);
//
// 				#ifdef DEBUG
// 					cout << "Orden actualizado " << _candidateOrders[iChannelOrder] << " " << processedParticle->GetChannelOrderAPP(iChannelOrder) << " ";
// 				#endif
// 			}

			#ifdef DEBUG
				cout << endl;
				cout << "El peso se actualiza multiplicando por " << sumLikelihoodsProd/proposal << endl;
			#endif

			// the weight is updated
			processedParticle->SetWeight((sumLikelihoodsProd/proposal)*processedParticle->GetWeight());

			#ifdef DEBUG2
				cout << "---------" << " " << iParticle << " " << "---------" << endl << Util::ToMatrix(sampledSmoothingVector,columnwise,_N) << "Peso actualizado por " << sumLikelihoodsProd/proposal << " | ";
				for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
				{
					cout << "Orden " << _candidateOrders[iChannelOrder] << ": " << processedParticle->GetChannelOrderAPP(iChannelOrder) << " | ";
				}
				cout << endl;
			#endif

			#ifdef DEBUG
				getchar();
			#endif
		} // for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)

		_particleFilter.NormalizeWeights();

		#ifdef DEBUG2
			tVector pesos = _particleFilter.GetWeightsVector();
			int iMaximoPeso;
			Util::Max(pesos,iMaximoPeso);
			cout << "La mejor partícula fue " << iMaximoPeso << " con peso " << pesos(iMaximoPeso) << endl;
		#endif

		#ifdef DEBUG2
			getchar();
		#endif

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
            _resamplingAlgorithm->Resample(&_particleFilter);
	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)

	delete[] newChannelOrderAPPs;
}

