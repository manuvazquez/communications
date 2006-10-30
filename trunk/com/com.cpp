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

#define HOSTNAME_LENGTH 50

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <unistd.h>
#include <string>

#include "types.h"
#include <Alphabet.h>
#include <Bits.h>
#include <ARprocess.h>
#include <ARchannel.h>
#include <SprawlingMemoryARMIMOChannel.h>
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
#include <UnknownChannelOrderAlgorithm.h>
#include <ML_UnknownChannelOrderSMCAlgorithm.h>
#include <ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm.h>
#include <ISIR.h>
#include <LinearFilterBasedISIRAlgorithm.h>
#include <LinearFilterBasedUnknownChannelOrderSMCAlgorithm.h>
#include <LinearFilterBasedUnknownChannelOrderModelSMCAlgorithm.h>

#include <ResamplingCriterion.h>
#include <StdResamplingAlgorithm.h>
#include <ByChannelOrderResamplingAlgorithm.h>
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

using namespace std;

int main(int argc,char* argv[])
{
    double pe,mse;

    // PARAMETERS
    int nFrames = 500;
    int L=3,N=2,m=2,K=300;
    int trainSeqLength = 30;
    int nParticles = 30;
    int d = m -1;
    char outputFileName[HOSTNAME_LENGTH+4] = "res_";

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

	// unknown channel order
	vector<int> candidateChannelOrders;
	candidateChannelOrders.push_back(2);candidateChannelOrders.push_back(3);candidateChannelOrders.push_back(4);/*candidateChannelOrders.push_back(5);*/

	// sprawling memory transition probability matrix
	tMatrix transitionProbabilitiesMatrix(candidateChannelOrders.size(),candidateChannelOrders.size());

	// probability of staying in the same channel order
	for(int i=0;i<candidateChannelOrders.size();i++)
		transitionProbabilitiesMatrix(i,i) = 0.7;
	transitionProbabilitiesMatrix(0,1) = 0.2;
	transitionProbabilitiesMatrix(0,2) = 0.1;
	transitionProbabilitiesMatrix(1,0) = 0.15;
	transitionProbabilitiesMatrix(1,2) = 0.15;
	transitionProbabilitiesMatrix(2,0) = 0.1;
	transitionProbabilitiesMatrix(2,1) = 0.2;

    // linear detectors parameters
    double forgettingFactorDetector = 0.98;

    // alphabet is defined
    vector<vector<tBit> > secuenciasBits(2,vector<tBit>(1));
    secuenciasBits[0][0] = 0; secuenciasBits[1][0] = 1;
    vector<tSymbol> simbolos(2);
    simbolos[0] = -1; simbolos[1] = 1;
    Alphabet pam2(1,2,secuenciasBits,simbolos);

	// -----------------------------------------------------------------------------

	// a vector that will contain the names of the algorithms
	vector<string> algorithmsNames;

	char hostname[HOSTNAME_LENGTH];
	gethostname(hostname,HOSTNAME_LENGTH);

	strcat(outputFileName,hostname);

    // a specific preamble is generated...
    tMatrix preambulo(N,10);
    preambulo = -1.0;

    tRange rAllSymbolRows(0,N-1);
    tRange rTrainingSequenceSymbolVectors(preambulo.cols(),preambulo.cols()+trainSeqLength-1);

    // the last "d" observations are used only for smoothing and don't result in any detected vector
    tRange rSymbolVectorsToComputePe(preambulo.cols()+trainSeqLength,K+preambulo.cols()-1);
//     tRange rSymbolVectorsToComputePe(100,K+preambulo.cols()-1);

    // channel estimators are constructed for the different algorithms
    KalmanEstimator kalmanEstimator(mediaInicial,ARcoefficients[0],ARvariance);
    RLSEstimator RLSestimator(mediaInicial,forgettingFactor);
    LMSEstimator LMSestimator(mediaInicial,muLMS);

	// vectors of channel estimators and linear detectors for unknown channel order algorithms
	vector<ChannelMatrixEstimator *> kalmanChannelEstimators;
	vector<RLSEstimator *> RLSchannelEstimators;
	vector<LinearDetector *> unknownChannelOrderLinearDetectors;
	for(int iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		kalmanChannelEstimators.push_back(new KalmanEstimator(LaGenMatDouble::zeros(L,N*candidateChannelOrders[iChannelOrder]),ARcoefficients[0],ARvariance));

		RLSchannelEstimators.push_back(new RLSEstimator(LaGenMatDouble::zeros(L,N*candidateChannelOrders[iChannelOrder]),forgettingFactor));

		unknownChannelOrderLinearDetectors.push_back(new RMMSEDetector(L*candidateChannelOrders[iChannelOrder],N*(2*candidateChannelOrders[iChannelOrder]-1),pam2.Variance(),forgettingFactorDetector,N*candidateChannelOrders[iChannelOrder]));
	}

    // linear filters
    RMMSEDetector RMMSEdetector(L*(d+1),N*(m+d),pam2.Variance(),forgettingFactorDetector,N*(d+1));

	// the maximum of the candidate channel orders is computed
	int maxCandidateOrder = candidateChannelOrders[0];
	for(int i=1;i<candidateChannelOrders.size();i++)
		if(candidateChannelOrders[i]>maxCandidateOrder)
			maxCandidateOrder = candidateChannelOrders[i];

	// the algorithms with the higher smoothing lag require
	int nSmoothingInstants = maxCandidateOrder-1;

	// the preamble that will be passed to the unknown channel order algorithms
	tMatrix unknownChannelOrderAlgorithmsPreamble(N,maxCandidateOrder-1);
	unknownChannelOrderAlgorithmsPreamble = -1.0;

    // always the same resampling criterion and algorithms
    ResamplingCriterion criterioRemuestreo(0.9);
    StdResamplingAlgorithm algoritmoRemuestreo(criterioRemuestreo);

    // resampling algorith for the case of unknown channel order
    ByChannelOrderResamplingAlgorithm unknownChannelOrderResamplingAlgorithm(criterioRemuestreo);

    // matrices for results
    tMatrix overallPeMatrix;
    tMatrix overallMseMatrix;

    // we don't want the same bits to be generated over and over
    Random bitsRandomGenerator(0);

    for(int iFrame=0;iFrame<nFrames;iFrame++)
    {
        // bits are generated ...
        Bits bitsTransmitir(N,K+nSmoothingInstants,bitsRandomGenerator);

        // ... and then modulated by means of the alphabet
        tMatrix simbolosTransmitir = Modulator::Modulate(bitsTransmitir,pam2);

        // the preamble is set before the symbols to be transmitted
        simbolosTransmitir = Util::Append(preambulo,simbolosTransmitir);

        tMatrix trainingSequence = simbolosTransmitir(rAllSymbolRows,rTrainingSequenceSymbolVectors);

        // an AR channel is generated
        ARchannel canal(N,L,m,simbolosTransmitir.cols(),channelMean,channelVariance,ARcoefficients,ARvariance);

		// a channel order varying AR channel is generated
// 		SprawlingMemoryARMIMOChannel canal(N,L,simbolosTransmitir.cols(),candidateChannelOrders,transitionProbabilitiesMatrix,0,channelMean,channelVariance,ARcoefficients,ARvariance);

		// noise is generated according to the channel
		ChannelDependentNoise ruido(&canal);

        for(int iSNR=0;iSNR<SNRs.size();iSNR++)
        {
            cout << "SNR = " << SNRs[iSNR] << " [Trama " << iFrame << "]..." << endl;

			// noise SNR is set
            ruido.SetSNR(SNRs[iSNR],pam2.Variance());

            // transmission
            tMatrix observaciones = canal.Transmit(simbolosTransmitir,ruido);

            // algorithms are created
            vector<Algorithm *> algorithms;

            // ----------------------- ALGORITHMS TO RUN ----------------------------

//             algorithms.push_back(new ML_SMCAlgorithm ("Detector suavizado optimo",pam2,L,N,K+preambulo.cols(),m,&kalmanEstimator,preambulo,d,nParticles,algoritmoRemuestreo));

//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("Filtro lineal LMS",pam2,L,N,K+preambulo.cols(),m,&LMSestimator,&RMMSEdetector,preambulo,d,nParticles,algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance));

//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("Filtro lineal RLS",pam2,L,N,K+preambulo.cols(),m,&RLSestimator,&RMMSEdetector,preambulo,d,nParticles,algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance));

            algorithms.push_back(new ViterbiAlgorithm("Viterbi",pam2,L,N,K+preambulo.cols(),canal,preambulo,d));

//             algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimator("Estimador de Kalman con simbolos conocidos",pam2,L,N,K+preambulo.cols(),m,&kalmanEstimator,preambulo,simbolosTransmitir));

//             algorithms.push_back(new ISIR("ISIR",pam2,L,N,K+preambulo.cols(),kalmanChannelEstimators,preambulo,preambulo.cols(),d,nParticles,&algoritmoRemuestreo));

//             algorithms.push_back(new LinearFilterBasedUnknownChannelOrderSMCAlgorithm("Linear Filter Unknown Channel Order",pam2,L,N,K+preambulo.cols(),kalmanChannelEstimators,unknownChannelOrderLinearDetectors,preambulo,preambulo.cols(),d,nParticles,&algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance,canal,simbolosTransmitir));

//             algorithms.push_back(new LinearFilterBasedISIRAlgorithm("Linear Filter Based ISIR with smoothing",pam2,L,N,K+preambulo.cols(),kalmanChannelEstimators,unknownChannelOrderLinearDetectors,preambulo,preambulo.cols(),d,nParticles,&algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance,canal,simbolosTransmitir));

//             algorithms.push_back(new LinearFilterBasedUnknownChannelOrderModelSMCAlgorithm("Linear Filter Based ISIR decidiendose por un modelo",pam2,L,N,K+preambulo.cols(),kalmanChannelEstimators,unknownChannelOrderLinearDetectors,preambulo,preambulo.cols(),d,nParticles,&algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance,canal,simbolosTransmitir));

			// -------------------- Los algoritmos de siempre considerando diferentes m's

			// m=2
//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("Filtro lineal RLS suponiendo m = 2",pam2,L,N,K+preambulo.cols(),candidateChannelOrders[0],RLSchannelEstimators[0],unknownChannelOrderLinearDetectors[0],preambulo,candidateChannelOrders[0]-1,nParticles,algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance));

			//m=3
//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("Filtro lineal RLS suponiendo m = 3",pam2,L,N,K+preambulo.cols(),candidateChannelOrders[1],RLSchannelEstimators[1],unknownChannelOrderLinearDetectors[1],preambulo,candidateChannelOrders[1]-1,nParticles,algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance));

			//m=4
// //             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("Filtro lineal RLS suponiendo m = 4",pam2,L,N,K+preambulo.cols(),candidateChannelOrders[2],RLSchannelEstimators[2],unknownChannelOrderLinearDetectors[2],preambulo,candidateChannelOrders[2]-1,nParticles,algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance));



// -----------------------------------------------------------------------------

//             algorithms.push_back(new ML_MultipleChannelEstimatorsPerParticleSMCAlgorithm("Estimador del orden del canal",pam2,L,N,K+preambulo.cols(),kalmanChannelEstimators,preambulo,preambulo.cols(),d,nParticles,&algoritmoRemuestreo));

// 			algorithms.push_back(new ML_UnknownChannelOrderSMCAlgorithm ("ML Unknown Channel Order",pam2,L,N,K+m-1,kalmanChannelEstimators,unknownChannelOrderAlgorithmsPreamble,m-1,d,nParticles,&unknownChannelOrderResamplingAlgorithm,&algoritmoRemuestreo,simbolosTransmitir));

//             algorithms.push_back(new ML_UnknownChannelOrderSMCAlgorithm ("ML Unknown Channel Order",pam2,L,N,K+unknownChannelOrderAlgorithmsPreamble.cols(),kalmanChannelEstimators,unknownChannelOrderAlgorithmsPreamble,m-1,d,nParticles,&unknownChannelOrderResamplingAlgorithm,&algoritmoRemuestreo,simbolosTransmitir));

            // ----------------------------------------------------------------------

            // here the number of algoriths is known. So, the first iteration:
            if(iFrame==0 && iSNR==0)
            {
                overallPeMatrix.resize(SNRs.size(),algorithms.size());
                overallPeMatrix = 0.0;

                overallMseMatrix.resize(SNRs.size(),algorithms.size());
                overallMseMatrix = 0.0;

				// we fill the vector with the names of the algorithms
				for(int iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
				{
					algorithmsNames.push_back(algorithms[iAlgorithm]->GetName());
				}
            }

            // algorithms are executed
            for(int iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
            {
                algorithms[iAlgorithm]->Run(observaciones,ruido.Variances(),trainingSequence);

                pe = algorithms[iAlgorithm]->SER(simbolosTransmitir(rAllSymbolRows,rSymbolVectorsToComputePe));
                mse = algorithms[iAlgorithm]->MSE(canal.Range(preambulo.cols(),K+preambulo.cols()-1));

                cout << algorithms[iAlgorithm]->GetName() << ": Pe = " << pe << " , MSE = " << mse << endl;

                // the error probability is accumulated
                overallPeMatrix(iSNR,iAlgorithm) += pe;

                // and the MSE
                overallMseMatrix(iSNR,iAlgorithm) += mse;

                delete algorithms[iAlgorithm];
            }
        } // for(int iSNR=0;iSNR<SNRs.size();iSNR++)


		// ----------------- VARIABLES SAVING ----------------------
		ofstream f(outputFileName,ofstream::trunc);

		tMatrix auxOverallPe = overallPeMatrix;
		auxOverallPe *= 1.0/(double)(iFrame+1);
		Util::MatrixToStream(auxOverallPe,"pe",f);

		tMatrix auxOverallMse = overallMseMatrix;
		auxOverallMse *= 1.0/(double)(iFrame+1);
		Util::MatrixToStream(auxOverallMse,"mse",f);

        Util::ScalarToStream(iFrame+1,"nFrames",f);

		Util::StringsVectorToStream(algorithmsNames,"algorithmsNames",f);
        Util::ScalarToStream(L,"L",f);
        Util::ScalarToStream(N,"N",f);
        Util::ScalarToStream(m,"m",f);
        Util::ScalarToStream(K,"K",f);
        Util::ScalarToStream(trainSeqLength,"trainSeqLength",f);
        Util::ScalarToStream(nParticles,"nParticles",f);
        Util::ScalarToStream(d,"d",f);
		Util::IntsVectorToStream(candidateChannelOrders,"candidateOrders",f);
		Util::IntsVectorToStream(SNRs,"SNRs",f);

		f.close();
		// ---------------------------------------------------------

    } // for(int iFrame=0;iFrame<nFrames;iFrame++)

    overallPeMatrix *= 1.0/nFrames;
    overallMseMatrix *= 1.0/nFrames;

    cout << "Overall SER:" << endl;
    Util::Print(overallPeMatrix);

    cout << "Overall MSE:" << endl;
    Util::Print(overallMseMatrix);
}
