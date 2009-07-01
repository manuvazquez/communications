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
#include "TriangularizationBasedSMCAlgorithm.h"

// #define DEBUG4
// #define DEBUG2

TriangularizationBasedSMCAlgorithm::TriangularizationBasedSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int frameLength, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances,double ARcoefficient,double ARprocessVariance): SMCAlgorithm(name, alphabet, L, N, frameLength, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances),_ARcoefficient(ARcoefficient),_ARprocessVariance(ARprocessVariance)
{
//     _randomParticlesInitilization = true;
}

void TriangularizationBasedSMCAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
    int iParticle,iSmoothing,iAlphabet,iSampled;
    double proposal,observationWithouNoise,sumProb,likelihoodsProd;
    vector<tMatrix> matricesToStack(_d+1,tMatrix(_L,_Nm));
    tRange rAllObservationsRows(0,_L-1),rAllSymbolRows(0,_N-1);
    tRange rAllStackedObservationsRows(0,_L*(_d+1)-1);
    tRange rFirstmMinus1symbolVectors(0,_m-2),rFirstmSymbolVectors(0,_m-1);
    tMatrix stackedChannelMatrixSubstract;
    tMatrix stackedChannelMatrixMinus(_L*(_d+1),_N*(_d+1)),stackedChannelMatrixMinusFlipped;
    tMatrix stackedChannelMatrixMinusFlippedTransposeStackedChannelMatrixMinusFlipped(_N*(_d+1),_N*(_d+1));
    tVector stackedObservationsMinus(_L*(_d+1));
    tMatrix L,invLstackedChannelMatrixMinusTrans(_N*(_d+1),_L*(_d+1));
    tMatrix observationsCovariance = LaGenMatDouble::zeros(_L*(_d+1),_L*(_d+1));
    tLongIntVector piv(_N*(_d+1));
    tVector transformedStackedObservationsMinus(_N*(_d+1));
    tMatrix transformedStackedObservationsCovariance(_N*(_d+1),_N*(_d+1));
    tMatrix invLstackedChannelMatrixMinusTransObservationsCovariance(_N*(_d+1),_L*(_d+1));
    tMatrix U(_N*(_d+1),_N*(_d+1));
    tMatrix involvedSymbolVectors = LaGenMatDouble::zeros(_N,_m+_d);
    int NmMinus1 = _N*(_m-1);
    tVector symbolProbabilities(_alphabet.length());

    for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
    {
#ifdef DEBUG
        cout << "iObservationToBeProcessed = " << iObservationToBeProcessed << endl;
#endif
        // already detected symbol vectors involved in the current detection
        tRange rAlreadyDetectedSymbolVectors(iObservationToBeProcessed-_m+1,iObservationToBeProcessed-1);

        // observation matrix columns that are involved in the smoothing
        tRange rSmoothingRange(iObservationToBeProcessed,iObservationToBeProcessed+_d);

        // the stacked observations vector
        tVector stackedObservations = Util::ToVector(observations(rAllObservationsRows,rSmoothingRange),columnwise);

        for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)
        {
            ParticleWithChannelEstimation *processedParticle = dynamic_cast <ParticleWithChannelEstimation *> (_particleFilter->GetParticle(iParticle));

            // the already detected symbol vectors are stored in "involvedSymbolVectors"
            involvedSymbolVectors(rAllSymbolRows,rFirstmMinus1symbolVectors).inject(processedParticle->GetSymbolVectors(rAlreadyDetectedSymbolVectors));

            // predicted channel matrices are stored in a vector in order to stack them
            // (first one is obtained via the Kalman Filter)
            matricesToStack[0] = (dynamic_cast<KalmanEstimator *> (processedParticle->GetChannelMatrixEstimator(_estimatorIndex)))->SampleFromPredictive();

            for(iSmoothing=1;iSmoothing<=_d;iSmoothing++)
            {
#ifdef DEBUG2
                cout << "iSmoothing = " << iSmoothing << endl;
#endif
                // matricesToStack[iSmoothing] = _ARcoefficient * matricesToStack[iSmoothing-1] + rand(_L,_Nm)*_ARprocessVariance
                Util::Add(matricesToStack[iSmoothing-1],StatUtil::RandnMatrix(_L,_Nm,0.0,_ARprocessVariance),matricesToStack[iSmoothing],_ARcoefficient,1.0);

                // "stackedChannelMatrixSubstract" will be used to substract the contribution of the already detected symbol vectors from the observations
//                 stackedChannelMatrixSubstract(tRange((iSmoothing-1)*_L,iSmoothing*_L-1),tRange(_N*(iSmoothing-1),(_m-1)*_N-1)).inject(matricesToStack[iSmoothing](rAllObservationsRows,tRange(0,(_m-iSmoothing)*_N-1)));
            }

            // matrices are stacked to give
            tMatrix stackedChannelMatrix = HsToStackedH(matricesToStack);

            stackedChannelMatrixMinus = stackedChannelMatrix(rAllStackedObservationsRows,tRange((_m-1)*_N,stackedChannelMatrix.cols()-1));

            stackedObservationsMinus = SubstractKnownSymbolsContribution(matricesToStack,_m,0,_d,stackedObservations,involvedSymbolVectors(rAllSymbolRows,rFirstmMinus1symbolVectors));

#ifdef DEBUG4
            cout << stackedObservationsMinus << endl;
            cout << rAllSymbolRows << endl << rFirstmMinus1symbolVectors << endl;
            cout << "con la funcion" << endl << SubstractKnownSymbolsContribution(matricesToStack,_m,0,_d,stackedObservations,involvedSymbolVectors(rAllSymbolRows,rFirstmMinus1symbolVectors)) << endl;
#endif

            // we want to start sampling the present symbol vector, not the future ones
            stackedChannelMatrixMinusFlipped = Util::flipLR(stackedChannelMatrixMinus);

#ifdef DEBUG2
            cout << "stackedChannelMatrixMinusFlipped" << endl << stackedChannelMatrixMinusFlipped;
#endif

            // stackedChannelMatrixMinusFlippedTransposeStackedChannelMatrixMinusFlipped = stackedChannelMatrixMinusFlipped'*stackedChannelMatrixMinusFlipped
            Blas_Mat_Trans_Mat_Mult(stackedChannelMatrixMinusFlipped,stackedChannelMatrixMinusFlipped,stackedChannelMatrixMinusFlippedTransposeStackedChannelMatrixMinusFlipped);

#ifdef DEBUG2
            cout << "Se calcula cholesky de" << endl << stackedChannelMatrixMinusFlippedTransposeStackedChannelMatrixMinusFlipped;
#endif

            // Cholesky decomposition is computed
            L = Util::cholesky(stackedChannelMatrixMinusFlippedTransposeStackedChannelMatrixMinusFlipped);

            // we also obtain the upper triangular matrix in U
            Util::transpose(L,U);

#ifdef DEBUG2
            cout << "la matriz de cholesky" << endl << L;
#endif

            // invL = inverse(L)
            tMatrix invL = L;
            LUFactorizeIP(invL,piv);
            LaLUInverseIP(invL,piv);

            // invLstackedChannelMatrixMinusTrans = invL*stackedChannelMatrixMinusFlipped'
            Blas_Mat_Mat_Trans_Mult(invL,stackedChannelMatrixMinusFlipped,invLstackedChannelMatrixMinusTrans);

            // the transformed observations are computed
            // transformedStackedObservationsMinus = invLstackedChannelMatrixMinusTrans * stackedObservationsMinus
            Blas_Mat_Vec_Mult(invLstackedChannelMatrixMinusTrans,stackedObservationsMinus,transformedStackedObservationsMinus);

            // the covariance of the transformed observations is computed...

            // ...starting by the covariance of the normal observations
            for(iSmoothing=0;iSmoothing<_d+1;iSmoothing++)
                for(int i=0;i<_L;i++)
                    observationsCovariance(iSmoothing*_L+i,iSmoothing*_L+i) = noiseVariances[iObservationToBeProcessed+iSmoothing];

            // invLstackedChannelMatrixMinusTransObservationsCovariance = invLstackedChannelMatrixMinusTrans * observationsCovariance
            Blas_Mat_Mat_Mult(invLstackedChannelMatrixMinusTrans,observationsCovariance,invLstackedChannelMatrixMinusTransObservationsCovariance);

            // transformedStackedObservationsCovariance = invLstackedChannelMatrixMinusTransObservationsCovariance * invLstackedChannelMatrixMinusTrans'
            Blas_Mat_Mat_Trans_Mult(invLstackedChannelMatrixMinusTransObservationsCovariance,invLstackedChannelMatrixMinusTrans,transformedStackedObservationsCovariance);

#ifdef DEBUG2
            cout << "transformedStackedObservationsCovariance" << endl << transformedStackedObservationsCovariance;
            cout << "Una tecla..."; getchar();
#endif

            // the evaluated proposal function (necessary for computing the weights) is initialized
            proposal = 1.0;

            for(int iSampledSymbol=_N*(_d+1)-1,iWithinMatrix=NmMinus1;iSampledSymbol>=0;iSampledSymbol--,iWithinMatrix++)
            {
                observationWithouNoise = 0.0;
                int jU,iS;
                for(jU=_N*(_d+1)-1,iS=NmMinus1;jU>iSampledSymbol;jU--,iS++)
                    observationWithouNoise += U(iSampledSymbol,jU)*involvedSymbolVectors(iS % _N,iS / _N);

                sumProb = 0.0;

                // the probability for each posible symbol alphabet is computed
                for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                {
#ifdef DEBUG2
                    cout << "jU = " << jU << " iSampledSymbol = " << iSampledSymbol << endl;
                    cout << "ehhh: " << observationWithouNoise+U(iSampledSymbol,jU)*_alphabet[iAlphabet] << endl;
#endif
                    symbolProbabilities(iAlphabet) = StatUtil::NormalPdf(transformedStackedObservationsMinus(iSampledSymbol),observationWithouNoise+U(iSampledSymbol,jU)*_alphabet[iAlphabet],transformedStackedObservationsCovariance(iSampledSymbol,iSampledSymbol));

                    // for normalization purposes
                    sumProb += symbolProbabilities(iAlphabet);
                }

                try {
                    for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                        symbolProbabilities(iAlphabet) /= sumProb;
                }catch(exception e){
                    cout << "TriangularizationBasedSMCAlgorithm::Process: the sum of the probabilities is null." << endl;
                    cout <<  __FILE__  << "(line " << __LINE__ << ") :" << endl;
                    for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                        symbolProbabilities(iAlphabet) = 1.0/double(_alphabet.length());
                }

#ifdef DEBUG
                cout << "Las probabilidades calculadas" << endl << symbolProbabilities;
#endif

                iSampled = StatUtil::discrete_rnd(symbolProbabilities);
                involvedSymbolVectors(iWithinMatrix % _N,iWithinMatrix / _N) = _alphabet[iSampled];
                proposal *= symbolProbabilities(iSampled);

            } // for(int iSampledSymbol=_N*(_d+1)-1,iWithinMatrix=NmMinus1;iSampledSymbol>=0;iSampledSymbol--,iWithinMatrix++)

            processedParticle->SetSymbolVector(iObservationToBeProcessed,involvedSymbolVectors.col(_m-1));

#ifdef DEBUG
            cout << "involvedSymbolVectors" << endl << involvedSymbolVectors;
#endif

            likelihoodsProd = Smoothedlikelihood(matricesToStack,involvedSymbolVectors,processedParticle,iObservationToBeProcessed,observations,noiseVariances);

            // the weight is updated
            processedParticle->SetWeight((likelihoodsProd/proposal)*processedParticle->GetWeight());

            // and the estimation of the channel matrix
            processedParticle->SetChannelMatrix(_estimatorIndex,iObservationToBeProcessed,
                                                processedParticle->GetChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),
                                                    involvedSymbolVectors(rAllSymbolRows,rFirstmSymbolVectors),noiseVariances[iObservationToBeProcessed]));

        } // for(iParticle=0;iParticle<_particleFilter->Capacity();iParticle++)

        _particleFilter->NormalizeWeights();

        // if it's not the last time instant
        if(iObservationToBeProcessed<(_K-1))
            _resamplingAlgorithm->ResampleWhenNecessary(_particleFilter);

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_K;iObservationToBeProcessed++)
}

