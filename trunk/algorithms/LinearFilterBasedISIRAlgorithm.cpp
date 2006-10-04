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
#include "LinearFilterBasedISIRAlgorithm.h"

LinearFilterBasedISIRAlgorithm::LinearFilterBasedISIRAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators,vector<LinearDetector *> linearDetectors, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_allObservationRows(0,_L-1),_linearDetectors(linearDetectors.size()),_particleFilter(nParticles),_ARcoefficient(ARcoefficient),_samplingVariance(samplingVariance),_ARprocessVariance(ARprocessVariance)
{
    if(linearDetectors.size()!=_candidateOrders.size())
        throw RuntimeException("LinearFilterBasedISIRAlgorithm::LinearFilterBasedISIRAlgorithm: nº of detectors and number of channel matrix estimators (and candidate orders) are different.");

    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _linearDetectors[iChannelOrder] = linearDetectors[iChannelOrder]->Clone();
}


LinearFilterBasedISIRAlgorithm::~LinearFilterBasedISIRAlgorithm()
{
    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        delete _linearDetectors[iChannelOrder];
}

vector<vector<tMatrix> > LinearFilterBasedISIRAlgorithm::ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
    int iChannelOrder;

    for(int i=_iFirstObservation;i<_iFirstObservation+trainingSequence.cols();i++)
    {
        for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        {
            // the observations from i to i+d are stacked
            tRange smoothingRange(i,i+_candidateOrders[iChannelOrder]-1);
            tVector stackedObservationsVector = Util::ToVector(observations(_allObservationRows,smoothingRange),columnwise);

            _linearDetectors[iChannelOrder]->StateStep(stackedObservationsVector);
        }
    }

    return UnknownChannelOrderAlgorithm::ProcessTrainingSequence(observations,noiseVariances,trainingSequence);
}

void LinearFilterBasedISIRAlgorithm::InitializeParticles()
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

void LinearFilterBasedISIRAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	int iParticle,iSmoothing,iRow,iSampledSymbol,iAlphabet,iSampled,iChannelOrder;
	int m,d,Nm;
	vector<vector<tMatrix> > matricesToStack(_candidateOrders.size());
	tVector sampledVector(_N),sampledSmoothingVector(_N*_maxOrder);
	double proposal,s2q,sumProb,likelihoodsProd,sumLikelihoodsProd,sumChannelOrderAPPs;

	// each matrix in symbolProb contains the probabilities connected to a channelOrder: symbolProb(i,j) is the p(i-th symbol=alphabet[j]). They are initialized with the matrix
	tMatrix symbolProbAux = LaGenMatDouble::zeros(_N*_maxOrder,_alphabet.Length());
// 	symbolProbAux *= 1.0/(double)_alphabet.Length(); // <--------------------------------------

	vector<tMatrix> symbolProb(_candidateOrders.size(),symbolProbAux);

	// "overallSymbolProb" will combine the previous probabilities accordin to the APP of the channel order
	tMatrix overallSymbolProb(_N*_maxOrder,_alphabet.Length());

	tMatrix forWeightUpdateNeededSymbols(_N,2*_maxOrder-1);
	tMatrix noiseCovariances[_maxOrder];
	tVector predictedNoiselessObservation(_L);

	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
// 		cout << "iObservationToBeProcessed: " << iObservationToBeProcessed << endl;

		// observation matrix columns that are involved in the smoothing
		tRange smoothingRange(iObservationToBeProcessed,iObservationToBeProcessed+_maxOrder-1);

		// the stacked observations vector
		tVector stackedObservations = Util::ToVector(observations(_allObservationRows,smoothingRange),columnwise);

		// stacked noise covariance needs to be constructed
		tMatrix stackedNoiseCovariance = LaGenMatDouble::zeros(_L*_maxOrder,_L*_maxOrder);
		for(iSmoothing=0;iSmoothing<_maxOrder;iSmoothing++)
			for(iRow=0;iRow<_L;iRow++)
				stackedNoiseCovariance(iSmoothing*_L+iRow,iSmoothing*_L+iRow) = noiseVariances[iObservationToBeProcessed+iSmoothing];

		// required noise covariances are computed from the noise variances
		for(iSmoothing=0;iSmoothing<=_maxOrder-1;iSmoothing++)
		{
			noiseCovariances[iSmoothing] = LaGenMatDouble::eye(_L);
			noiseCovariances[iSmoothing] *= noiseVariances[iObservationToBeProcessed+iSmoothing];
		}

		for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
		{
			ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP *>(_particleFilter.GetParticle(iParticle));

			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				m = _candidateOrders[iChannelOrder];
				d = m-1;
				Nm = _N*m;
				matricesToStack[iChannelOrder] = vector<tMatrix>(_maxOrder,tMatrix(_L,Nm));
				tRange rPresentChannelOrderObservations(0,_L*(d+1)-1);
				tRange rStackedNoiseCovariance(0,_L*(d+1)-1);
				tMatrix s2qAux(_L*(d+1),_L*(d+1));
				tVector s2qAuxFilter(_L*(d+1));

				// predicted channel matrices are stored in a vector in order to stack them

				// matricesToStack[iChannelOrder][0] = _ARcoefficient * <lastEstimatedChannelMatrix> + randn(_L,Nm)*_samplingVariance
				Util::Add((processedParticle->GetChannelMatrixEstimator(iChannelOrder))->LastEstimatedChannelMatrix(),StatUtil::RandnMatrix(_L,Nm,0.0,_samplingVariance),matricesToStack[iChannelOrder][0],_ARcoefficient,1.0);

				for(iSmoothing=1;iSmoothing<_maxOrder;iSmoothing++)
				{
					// matricesToStack[iChannelOrder][iSmoothing] = _ARcoefficient * matricesToStack[iChannelOrder][iSmoothing-1] + rand(_L,Nm)*_ARprocessVariance
					Util::Add(matricesToStack[iChannelOrder][iSmoothing-1],StatUtil::RandnMatrix(_L,Nm,0.0,_ARprocessVariance),matricesToStack[iChannelOrder][iSmoothing],_ARcoefficient,1.0);
				}

				// matrices are stacked to give
				tMatrix stackedChannelMatrix = HsToStackedH(matricesToStack[iChannelOrder],m,d);

				// the estimated stacked channel matrix is used to obtain soft estimations
				// of the transmitted symbols
				tVector softEstimations = processedParticle->GetLinearDetector(iChannelOrder)->Detect(stackedObservations(rPresentChannelOrderObservations),stackedChannelMatrix);

				tMatrix filter = processedParticle->GetLinearDetector(iChannelOrder)->ComputedFilter();

				// operations needed to computed the sampling variance

				//s2qAux = _alphabet.Variance() * stackedChannelMatrix * stackedChannelMatrix^H
				Blas_Mat_Mat_Trans_Mult(stackedChannelMatrix,stackedChannelMatrix,s2qAux,_alphabet.Variance());

				// s2qAux = s2qAux + stackedNoiseCovariance
				Util::Add(s2qAux,stackedNoiseCovariance(rStackedNoiseCovariance,rStackedNoiseCovariance),s2qAux);

// 				cout << "Antes del calculo de probabilidades" << endl;

				// sampling
				for(iSampledSymbol=0;iSampledSymbol<(_N*(d+1));iSampledSymbol++)
				{
					// s2qAuxFilter = s2qAux * filter.col(iSampledSymbol)
					Blas_Mat_Vec_Mult(s2qAux,filter.col(iSampledSymbol),s2qAuxFilter);


					s2q = _alphabet.Variance()*(1.0 - 2.0*Blas_Dot_Prod(filter.col(iSampledSymbol),stackedChannelMatrix.col(iSampledSymbol))) + Blas_Dot_Prod(filter.col(iSampledSymbol),s2qAuxFilter);

					double sumProb = 0.0;
					// the probability for each posible symbol alphabet is computed
					for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
					{
						symbolProb[iChannelOrder](iSampledSymbol,iAlphabet) = StatUtil::NormalPdf(softEstimations(iSampledSymbol),_alphabet[iAlphabet],s2q);

						// the computed pdf is accumulated for normalizing purposes
						sumProb += symbolProb[iChannelOrder](iSampledSymbol,iAlphabet);
					}

					try {
						for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
							symbolProb[iChannelOrder](iSampledSymbol,iAlphabet) /= sumProb;
					}catch(exception e){
						cout << "The sum of the probabilities is null." << endl;
						for(iAlphabet=0;iAlphabet<_alphabet.Length();iAlphabet++)
							symbolProb[iChannelOrder](iSampledSymbol,iAlphabet) = 0.5;
					}
				}
// 				cout << "Despues del calculo de probabilidades" << endl;

// 				cout << "Matriz de probabilidades para este orden" << endl << symbolProb[iChannelOrder] << endl;
			}

// 			for(iChannelOrder=1;iChannelOrder<_candidateOrders.size();iChannelOrder++)
// 			{
// 				cout << "Para iChannelOrder " << iChannelOrder << " " << processedParticle->GetChannelOrderAPP(iChannelOrder);
// 			}

			// the probabilities of the different channel orders are weighted according to the a posteriori probability of the channel order in the previous time instant
			overallSymbolProb = symbolProb[0];
			overallSymbolProb *= processedParticle->GetChannelOrderAPP(0);
			for(iChannelOrder=1;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				Util::Add(overallSymbolProb,symbolProb[iChannelOrder],overallSymbolProb,1.0,processedParticle->GetChannelOrderAPP(iChannelOrder));
			}

			proposal = 1.0;

// 			cout << "Antes de samplear" << endl;

			// the symbols are sampled from the above combined probabilities
			for(iSampledSymbol=0;iSampledSymbol<_N*_maxOrder;iSampledSymbol++)
			{
				int iSampled = (StatUtil::Discrete_rnd(1,overallSymbolProb.row(iSampledSymbol)))[0];
				sampledSmoothingVector(iSampledSymbol) = _alphabet[iSampled];

				proposal *= overallSymbolProb(iSampledSymbol,iSampled);
			}

// 			cout << "Antes de samplear" << endl;

			// sampled symbol vector is stored for the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledSmoothingVector(_allSymbolsRows));

// 			cout << "Guardados los símbolos" << endl;

			sumLikelihoodsProd = 0.0;
			sumChannelOrderAPPs = 0.0;

			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				m = _candidateOrders[iChannelOrder];
				d = m-1;

// 				cout << "Antes de lo de los rangos" << endl;

				// already detected symbol vectors involved in the current detection
				tRange rAlreadyDetectedSymbolVectors(iObservationToBeProcessed-m+1,iObservationToBeProcessed-1);

				tRange r0mMinus2(0,m-2),rSampledSymbolVectors(m-1,m+_maxOrder-2);
				tRange rFirstmSymbolVectors(0,m-1);

				// all the symbol vectors involved in the smoothing are kept in "forWeightUpdateNeededSymbols"
				// i) the already known:
				forWeightUpdateNeededSymbols(_allSymbolsRows,r0mMinus2).inject(processedParticle->GetSymbolVectors(rAlreadyDetectedSymbolVectors));

// 				cout << "Lo que muestreó" << endl << Util::ToMatrix(sampledSmoothingVector,columnwise,_N) << endl;

// 				cout << "Y lo va a meter en" << endl << forWeightUpdateNeededSymbols << endl;

// 				cout << "Pasado el primero" << endl;


				// ii) the just sampled
				forWeightUpdateNeededSymbols(_allSymbolsRows,rSampledSymbolVectors).inject(Util::ToMatrix(sampledSmoothingVector,columnwise,_N));

// 				cout << "Antes de empezar el calculo de la versosimilitudes" << endl;

// 				likelihoodsProd = 1.0;
				likelihoodsProd = processedParticle->GetChannelOrderAPP(iChannelOrder);

				for(iSmoothing=0;iSmoothing<=_maxOrder-1;iSmoothing++)
				{
					tRange rSymbolVectors(iSmoothing,iSmoothing+m-1);
					tVector stackedSymbolVector = Util::ToVector(forWeightUpdateNeededSymbols(_allSymbolsRows,rSymbolVectors),columnwise);

					// predictedNoiselessObservation = matricesToStack[iChannelOrder][iSmoothing] * stackedSymbolVector
					Blas_Mat_Vec_Mult(matricesToStack[iChannelOrder][iSmoothing],stackedSymbolVector,predictedNoiselessObservation);

					likelihoodsProd *= StatUtil::NormalPdf(observations.col(iObservationToBeProcessed+iSmoothing),predictedNoiselessObservation,noiseCovariances[iSmoothing]);

					// unnormalized APP channel order
					if(iSmoothing==0)
					{
						processedParticle->SetChannelOrderAPP(likelihoodsProd,iChannelOrder);

						// it's accumulated for normalization purposes
						sumChannelOrderAPPs += likelihoodsProd;
					}
				}

				// the estimation of the channel matrix is updated
				processedParticle->SetChannelMatrix(iChannelOrder,iObservationToBeProcessed, (processedParticle->GetChannelMatrixEstimator(iChannelOrder))->NextMatrix(observations.col(iObservationToBeProcessed),forWeightUpdateNeededSymbols(_allSymbolsRows,rFirstmSymbolVectors),noiseVariances[iObservationToBeProcessed]));
			}

			// the computed likelihood is accumulated
			sumLikelihoodsProd += likelihoodsProd;

			// the channel order APPs are normalized for the next iteration
			if(sumChannelOrderAPPs==0)
				throw RuntimeException("LinearFilterBasedISIRAlgorithm::Process: the sum of the channel order APP is zero.");
			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				processedParticle->SetChannelOrderAPP(processedParticle->GetChannelOrderAPP(iChannelOrder)/sumChannelOrderAPPs,iChannelOrder);
			}

			// the weight is updated
			processedParticle->SetWeight((likelihoodsProd/proposal)*processedParticle->GetWeight());
		}
// 		cout << "Antes de normalizar" << endl;
		_particleFilter.NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
            _resamplingAlgorithm->Resample(&_particleFilter);
	}
}

