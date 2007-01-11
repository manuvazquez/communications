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
#include "CME_UCO_SIS.h"

#define DEBUG

CME_UCO_SIS::CME_UCO_SIS(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, vector< LinearDetector * > linearDetectors, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, double ARcoefficient, double samplingVariance, double ARprocessVariance): UCO_SIS(name, alphabet, L, N, K, channelEstimators, linearDetectors, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm, ARcoefficient, samplingVariance, ARprocessVariance)
{
}


CME_UCO_SIS::~CME_UCO_SIS()
{
}


void CME_UCO_SIS::Process(const tMatrix& observations, vector< double > noiseVariances)
{
	int iParticle,iSmoothing,iRow,iSampledSymbol,iAlphabet,iSampled,iChannelOrder;
	int m,d,Nm,nLinearFiltersNeeded,iLinearFilterNeeded;
	vector<vector<tMatrix> > matricesToStack(_candidateOrders.size());
	tVector sampledVector(_N),sampledSmoothingVector(_N*_maxOrder);
	double proposal,s2q,sumProb,likelihoodsProd,sumLikelihoodsProd,channelOrderAPPsNormConstant;
	double *newChannelOrderAPPs = new double[_candidateOrders.size()];

	// each matrix in "symbolProb" contains the probabilities connected to a channelOrder: symbolProb(i,j) is the p(i-th symbol=alphabet[j]). They are initialized with zeros
	vector<tMatrix> symbolProb(_candidateOrders.size(),LaGenMatDouble::zeros(_N*_maxOrder,_alphabet.Length()));

	// "overallSymbolProb" will combine the previous probabilities accordin to the APP of the channel order
	tMatrix overallSymbolProb(_N*_maxOrder,_alphabet.Length());

	// 2*_maxOrder-1 = m_{max} + d_{max}
	tMatrix forWeightUpdateNeededSymbols(_N,2*_maxOrder-1);

	// _maxOrder = d_{max} + 1
	tMatrix noiseCovariances[_maxOrder];

	tVector predictedNoiselessObservation(_L);

	// ------------------ CME ADDON -------------------------

// 	vector<int> CMEstackedInstantsNeeded(_candidateOrders.size());
// 	vector<int> nTimeInstantsToStack(_candidateOrders.size());
//
// 	for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
// 	{
//
// 		// we want to know how many observations vectors we need to stack for applying CME
// 		CMEstackedInstantsNeeded[iChannelOrder] = (_N*_candidateOrders[iChannelOrder] - _L)/(_L - _N) + 1;
//
// 		// at least, we are going to stack "_maxOrder-1"
// 		nTimeInstantsToStack[iChannelOrder] = CMEstackedInstantsNeeded[iChannelOrder] > (_maxOrder-1) ? CMEstackedInstantsNeeded[iChannelOrder]:(_maxOrder-1);
//
// 		#ifdef DEBUG
// 			cout << "orden: " << _candidateOrders[iChannelOrder] << " requeridos por CME: " << CMEstackedInstantsNeeded[iChannelOrder] << " entonces: " << nTimeInstantsToStack[iChannelOrder] << endl;
// 		#endif
// 	}

	// we want to know how many observations vectors we need to stack for applying CME
	int nTimeInstantsToStackCME = (_N*_maxOrder - _L)/(_L - _N) + 1;
	nTimeInstantsToStackCME = nTimeInstantsToStackCME > 0 ? nTimeInstantsToStackCME:0;

	// at least, we are going to stack "_maxOrder-1"
	int nTimeInstantsToStack = nTimeInstantsToStackCME > (_maxOrder-1) ? nTimeInstantsToStackCME:(_maxOrder-1);

	#ifdef DEBUG
		cout << "nTimeInstantsToStackCME: " << nTimeInstantsToStackCME << " nTimeInstantsToStack: " << nTimeInstantsToStack << endl;
	#endif

	// ------------------------------------------------------

	for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
	{
		// observation matrix columns that are involved
		tRange rSmoothingRange(iObservationToBeProcessed,iObservationToBeProcessed+nTimeInstantsToStack);

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

		for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)
		{
			ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderAPP *>(_particleFilter.GetParticle(iParticle));

			for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			{
				m = _candidateOrders[iChannelOrder];
				d = m-1;
				Nm = _N*m;
				matricesToStack[iChannelOrder] = vector<tMatrix>(nTimeInstantsToStack+1,tMatrix(_L,Nm));

				tMatrix s2qAux(_L*(d+1),_L*(d+1));
				tVector s2qAuxFilter(_L*(d+1));

				// predicted channel matrices are sampled and stored in a vector in order to stack them

				// matricesToStack[iChannelOrder][0] = _ARcoefficient * <lastEstimatedChannelMatrix> + randn(_L,Nm)*_samplingVariance
				Util::Add(processedParticle->GetChannelMatrixEstimator(iChannelOrder)->LastEstimatedChannelMatrix(),StatUtil::RandnMatrix(_L,Nm,0.0,_samplingVariance),matricesToStack[iChannelOrder][0],_ARcoefficient,1.0);

				for(iSmoothing=1;iSmoothing<=nTimeInstantsToStack;iSmoothing++)
				{
					// matricesToStack[iChannelOrder][iSmoothing] = _ARcoefficient * matricesToStack[iChannelOrder][iSmoothing-1] + rand(_L,Nm)*_ARprocessVariance
					Util::Add(matricesToStack[iChannelOrder][iSmoothing-1],StatUtil::RandnMatrix(_L,Nm,0.0,_ARprocessVariance),matricesToStack[iChannelOrder][iSmoothing],_ARcoefficient,1.0);
				}

				// ------------------------------ CME ---------------------------------

				tMatrix CMEstackedChannelMatrix = HsToStackedH(matricesToStack[iChannelOrder],_candidateOrders[iChannelOrder],nTimeInstantsToStackCME);

				#ifdef DEBUG
					cout << "Orden de canal: " << _candidateOrders[iChannelOrder] << endl;
					cout << "La matriz apilada para CME es" << endl << CMEstackedChannelMatrix << endl;
					cout << "filas: " << CMEstackedChannelMatrix.rows() << " columnas: " << CMEstackedChannelMatrix.cols() << endl;
					cout << "debería tener " << _N*(_candidateOrders[iChannelOrder]+nTimeInstantsToStackCME) << " columnas" << endl;
				#endif

				int nRowsCMEstackedChannelMatrix = CMEstackedChannelMatrix.rows();
				int nColsCMEstackedChannelMatrix = CMEstackedChannelMatrix.cols();

				tMatrix HTH(nColsCMEstackedChannelMatrix,nColsCMEstackedChannelMatrix);
				// HTH = matricesToStack[iChannelOrder][0]^T * matricesToStack[iChannelOrder][0]
				Blas_Mat_Trans_Mat_Mult(CMEstackedChannelMatrix,CMEstackedChannelMatrix,HTH);

				// both the determinant and the inverse of HTH are computed
				tLongIntVector piv(nColsCMEstackedChannelMatrix);
				LUFactorizeIP(HTH,piv);

				double detHTH = 1.0;
				for(int i=0;i<nRowsCMEstackedChannelMatrix;i++)
					detHTH *= HTH(i,i);

				// because HTH is positive definite, the determinant is known to be positive
				detHTH = fabs(detHTH);

				#ifdef DEBUG
					cout << "El determinante sin denominador es " << detHTH << endl;
				#endif

				// inverse is computed
				LaLUInverseIP(HTH,piv);

				#ifdef DEBUG
					cout << "La inversa es " << endl << HTH << endl;
				#endif

				tMatrix HinvHTH(nRowsCMEstackedChannelMatrix,nColsCMEstackedChannelMatrix);
				// HinvHTH = CMEstackedChannelMatrix * inv(HTH)
				Blas_Mat_Mat_Mult(CMEstackedChannelMatrix,HTH,HinvHTH);

				tMatrix identityMinusHinvHTH_HT(nRowsCMEstackedChannelMatrix,nRowsCMEstackedChannelMatrix);
				// identityMinusHinvHTH_HT = HinvHTH * CMEstackedChannelMatrix
				// next, identity matrix will be added
				Blas_Mat_Mat_Trans_Mult(HinvHTH,CMEstackedChannelMatrix,identityMinusHinvHTH_HT,-1.0);

				// identityMinusHinvHTH_HT = I + identityMinusHinvHTH_HT
				for(int i=0;i<nRowsCMEstackedChannelMatrix;i++)
					identityMinusHinvHTH_HT(i,i)+= 1.0;

				tVector xTidentityMinusHinvHTH_HT(nRowsCMEstackedChannelMatrix);
				// xTidentityMinusHinvHTH_HT = stackedObservations(tRange(0,nRowsCMEstackedChannelMatrix-1))^T * identityMinusHinvHTH_HT = identityMinusHinvHTH_HT^T * stackedObservations(tRange(0,nRowsCMEstackedChannelMatrix-1))
				Blas_Mat_Trans_Vec_Mult(identityMinusHinvHTH_HT,stackedObservations(tRange(0,nRowsCMEstackedChannelMatrix-1)),xTidentityMinusHinvHTH_HT);

				double estimatedNoiseVariance = Blas_Dot_Prod(xTidentityMinusHinvHTH_HT,stackedObservations(tRange(0,nRowsCMEstackedChannelMatrix-1)));

				#ifdef DEBUG
					cout << "varianza estimada: " << estimatedNoiseVariance << endl;
					getchar();
				#endif

				double CME = (double)((nRowsCMEstackedChannelMatrix - nColsCMEstackedChannelMatrix)/2)*estimatedNoiseVariance/noiseVariances[iObservationToBeProcessed] + 0.5*log(detHTH/pow(2*M_PI*noiseVariances[iObservationToBeProcessed],nRowsCMEstackedChannelMatrix));

				#ifdef DEBUG
					cout << " CME: " << CME << endl;
					getchar();
				#endif

				// --------------------------------------------------------------------

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

				// the clone is dismissed
				delete linearDetectorBeingProccessed;

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
			for(iSampledSymbol=0;iSampledSymbol<_N*_maxOrder;iSampledSymbol++)
			{
				int iSampled = StatUtil::Discrete_rnd(overallSymbolProb.row(iSampledSymbol));

				sampledSmoothingVector(iSampledSymbol) = _alphabet[iSampled];

				proposal *= overallSymbolProb(iSampledSymbol,iSampled);
			}

			// sampled symbol vector is stored for the corresponding particle
			processedParticle->SetSymbolVector(iObservationToBeProcessed,sampledSmoothingVector(_allSymbolsRows));

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

				// all the symbol vectors involved in the smoothing are kept in "forWeightUpdateNeededSymbols":
				// i) the already known:
				forWeightUpdateNeededSymbols(_allSymbolsRows,r0mMinus2).inject(processedParticle->GetSymbolVectors(rAlreadyDetectedSymbolVectors));

				// ii) the just sampled
				forWeightUpdateNeededSymbols(_allSymbolsRows,rSampledSymbolVectors).inject(Util::ToMatrix(sampledSmoothingVector,columnwise,_N));

				likelihoodsProd = processedParticle->GetChannelOrderAPP(iChannelOrder);

				for(iSmoothing=0;iSmoothing<_maxOrder;iSmoothing++)
				{
					tRange rSymbolVectors(iSmoothing,iSmoothing+m-1);
					tVector stackedSymbolVector = Util::ToVector(forWeightUpdateNeededSymbols(_allSymbolsRows,rSymbolVectors),columnwise);

					// predictedNoiselessObservation = matricesToStack[iChannelOrder][iSmoothing] * stackedSymbolVector
					Blas_Mat_Vec_Mult(matricesToStack[iChannelOrder][iSmoothing],stackedSymbolVector,predictedNoiselessObservation);

					likelihoodsProd *= StatUtil::NormalPdf(observations.col(iObservationToBeProcessed+iSmoothing),predictedNoiselessObservation,noiseCovariances[iSmoothing]);

					// unnormalized APP channel order
					if(iSmoothing==0)
					{
						newChannelOrderAPPs[iChannelOrder] = likelihoodsProd;

						// it's accumulated for normalization purposes
						channelOrderAPPsNormConstant += likelihoodsProd;
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

			// the weight is updated
			processedParticle->SetWeight((sumLikelihoodsProd/proposal)*processedParticle->GetWeight());

		} // for(iParticle=0;iParticle<_particleFilter.Nparticles();iParticle++)

		_particleFilter.NormalizeWeights();

		// if it's not the last time instant
		if(iObservationToBeProcessed<(_K-1))
            _resamplingAlgorithm->Resample(&_particleFilter);
	} // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)

	delete[] newChannelOrderAPPs;
}

