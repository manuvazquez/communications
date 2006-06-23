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

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include "types.h"
#include <Alphabet.h>
#include <Bits.h>
#include <ARprocess.h>
#include <ARchannel.h>
#include <ChannelDependentNoise.h>
#include <Modulator.h>
#include <Demodulator.h>
#include <KalmanFilter.h>
#include <KalmanEstimator.h>
#include <RLSEstimator.h>
#include <LMSEstimator.h>
#include <RMMSEDetector.h>
#include <ML_SMCAlgorithm.h>
#include <LinearFilterBasedSMCAlgorithm.h>
#include <ViterbiAlgorithm.h>
#include <KnownSymbolsKalmanBasedChannelEstimator.h>
#include <ResamplingCriterion.h>
#include <StdResamplingAlgorithm.h>
#include <StatUtil.h>
#include <Util.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>
#include <Particle.h>
#include <ParticleWithChannelEstimation.h>


// #define MAX_NUMBER_ALGORITHMS 50
// #define EXPECTED_NUMBER_FRAMES 20

using namespace std;

int main(int argc,char* argv[])
{
    double pe,mse;

    // PARAMETERS
    int nFrames = 500;
    int L=3,N=2,m=2,K=300;
    int longSecEntr = 30;
    int nParticles = 30;
    int d = m -1;

    // SNRs to be processed
    vector<int> SNRs;
    SNRs.push_back(3);SNRs.push_back(6);SNRs.push_back(9);SNRs.push_back(12);SNRs.push_back(15);

    // AR process parameters
    vector<double> ARcoefficients(1);
    ARcoefficients[0] = 0.99999;
    double ARvariance=0.0001;

    // channel parameters
    double channelMean=0.0,channelVariance=1.0;

    // channel estimator parameters
    double samplingVariance = 0.0625;
    tMatrix mediaInicial(L,N*m);
    mediaInicial = 0.0;
    double forgettingFactor = 0.99;
    double muLMS = 0.05;

    // linear detectors parameters
    double forgettingFactorDetector = 0.98;

    // alphabet is defined
    vector<vector<tBit> > secuenciasBits(2,vector<tBit>(1));
    secuenciasBits[0][0] = 0; secuenciasBits[1][0] = 1;
    vector<tSymbol> simbolos(2);
    simbolos[0] = -1; simbolos[1] = 1;
    Alphabet pam2(1,2,secuenciasBits,simbolos);


    tRange rAllSymbolRows(0,N-1);
    tRange rTrainingSequenceSymbolVectors(m-1,m+longSecEntr-2);

    // channel estimators are constructed for the different algorithms
    KalmanEstimator kalmanEstimator(mediaInicial,ARcoefficients[0],ARvariance);
    RLSEstimator RLSestimator(mediaInicial,forgettingFactor);
    LMSEstimator LMSestimator(mediaInicial,muLMS);

    // linear filters
    RMMSEDetector RMMSEdetector(L*(d+1),N*(m+d),pam2.Variance(),forgettingFactorDetector,N*(d+1));

    // always the same resampling criterion and algorithms
    ResamplingCriterion criterioRemuestreo(0.9);
    StdResamplingAlgorithm algoritmoRemuestreo;

    // matrices for results
    tMatrix overallPeMatrix;
    tMatrix overallMseMatrix;

    for(int iFrame=0;iFrame<nFrames;iFrame++)
    {
        // bits are generated ...
        Bits bitsTransmitir(N,K);
    
        // ... and then modulated by means of the alphabet
        tMatrix simbolosTransmitir = Modulator::Modulate(bitsTransmitir,pam2);
    
        // a specific preamble is generated...
        tMatrix preambulo(N,m-1);
        preambulo = -1.0;
    
        // ...and set before the symbols to be transmitted
        simbolosTransmitir = Util::Append(preambulo,simbolosTransmitir);

        // notice that the last "d" observations are used only for smoothing and don't result in any detected vector
        tRange rSymbolVectorsToComputePe(m-1+longSecEntr,simbolosTransmitir.cols()-d-1);

        tMatrix trainingSequence = simbolosTransmitir(rAllSymbolRows,rTrainingSequenceSymbolVectors);

        // an AR channel is generated
        ARchannel canal(N,L,m,simbolosTransmitir.cols(),channelMean,channelVariance,ARcoefficients,ARvariance);

        for(int iSNR=0;iSNR<SNRs.size();iSNR++)
        {
            cout << "SNR = " << SNRs[iSNR] << " [Trama " << iFrame << "]..." << endl;
            
            // noise is generated
            ChannelDependentNoise ruido(&canal);
            ruido.SetSNR(SNRs[iSNR],pam2.Variance());

            // transmission
            tMatrix observaciones = canal.Transmit(simbolosTransmitir,ruido);

            // algorithms are created
            vector<Algorithm *> algorithms;

            // ----------------------- ALGORITHMS TO RUN ----------------------------

            algorithms.push_back(new ML_SMCAlgorithm ("Detector suavizado optimo",pam2,&kalmanEstimator,preambulo,m-1,nParticles,criterioRemuestreo,algoritmoRemuestreo));

            algorithms.push_back(new LinearFilterBasedSMCAlgorithm("Filtro lineal LMS",pam2,&LMSestimator,&RMMSEdetector,preambulo,m-1,nParticles,criterioRemuestreo,algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance));

            algorithms.push_back(new LinearFilterBasedSMCAlgorithm("Filtro lineal RLS",pam2,&RLSestimator,&RMMSEdetector,preambulo,m-1,nParticles,criterioRemuestreo,algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance));

            algorithms.push_back(new ViterbiAlgorithm("Viterbi",pam2,canal,preambulo,d));

            // we don't want the channel matrices corresponding to the smoothing observations to be detected, so we don't pass all the transmitted symbol vectors to the constructor
            algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimator("Estimador de Kalman con simbolos conocidos",pam2,&kalmanEstimator,preambulo,simbolosTransmitir(rAllSymbolRows,tRange(0,simbolosTransmitir.cols()-d-1))));

            // ----------------------------------------------------------------------

            KalmanEstimator kalmanEstimatorCopy(kalmanEstimator);

            // here the number of algoriths is known. So, the first iteration:
            if(iFrame==0 && iSNR==0)
            {
                overallPeMatrix.resize(SNRs.size(),algorithms.size());
                overallPeMatrix = 0.0;

                overallMseMatrix.resize(SNRs.size(),algorithms.size());
                overallMseMatrix = 0.0;
            }

            // now executed
            for(int iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
            {
                algorithms[iAlgorithm]->Run(observaciones,ruido.Variances(),trainingSequence);

                pe = algorithms[iAlgorithm]->SER(simbolosTransmitir(rAllSymbolRows,rSymbolVectorsToComputePe));
                mse = algorithms[iAlgorithm]->MSE(canal.Range(m-1+longSecEntr,simbolosTransmitir.cols()-d-1));

                cout << algorithms[iAlgorithm]->GetName() << ": Pe = " << pe << " , MSE = " << mse << endl;

                // the error probability is accumulated
                overallPeMatrix(iSNR,iAlgorithm) += pe;

                // and the MSE
                overallMseMatrix(iSNR,iAlgorithm) += mse;

                delete algorithms[iAlgorithm];
            }
        }
		ofstream f("res",ofstream::trunc);

		tMatrix auxOverallPe = overallPeMatrix;
		auxOverallPe *= 1.0/(double)(iFrame+1);
		Util::MatrixToStream(auxOverallPe,"pe",f);

		tMatrix auxOverallMse = overallMseMatrix;
		auxOverallMse *= 1.0/(double)(iFrame+1);
		Util::MatrixToStream(auxOverallMse,"mse",f);

		f.close();
    }

    overallPeMatrix *= 1.0/nFrames;
    overallMseMatrix *= 1.0/nFrames;

    cout << "Overall SER:" << endl;
    Util::Print(overallPeMatrix);

    cout << "Overall MSE:" << endl;
    Util::Print(overallMseMatrix);
}
