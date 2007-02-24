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

// #define DEBUG

#define HOSTNAME_LENGTH 50
#define DATE_LENGTH 100
#define SPRINTF_BUFFER 30

// the seed used to create the random objects is generated from the system time
// #define RANDOM_SEED

// wether or not, data regarding the channel orders APP evolution is saved
#define CHANNELORDERSAPP_SAVING

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <unistd.h>
#include <string>
#include <sys/time.h>

#include "types.h"
#include <Alphabet.h>
#include <Bits.h>
#include <ARprocess.h>

#include <ARchannel.h>
#include <ARoneChannelOrderPerTransmitAtennaMIMOChannel.h>
#include <SparklingMemoryARMIMOChannel.h>

#include <ChannelDependentNoise.h>
#include <NullNoise.h>
#include <Modulator.h>
#include <Demodulator.h>
#include <KalmanFilter.h>

#include <KalmanEstimator.h>
#include <RLSEstimator.h>
#include <LMSEstimator.h>
#include <OneChannelOrderPerTransmitAtennaWrapperEstimator.h>
#include <APPbasedChannelOrderEstimator.h>

#include <RMMSEDetector.h>
#include <MMSEDetector.h>

#include <DSISoptAlgorithm.h>
#include <LinearFilterBasedSMCAlgorithm.h>
#include <ViterbiAlgorithm.h>
#include <KnownSymbolsKalmanBasedChannelEstimator.h>
#include <UnknownChannelOrderAlgorithm.h>
#include <ISIS.h>
#include <USIS.h>

#include <Particle.h>
#include <ParticleWithChannelEstimation.h>
#include <ResamplingCriterion.h>
#include <StdResamplingAlgorithm.h>
#include <ByChannelOrderResamplingAlgorithm.h>
#include <StatUtil.h>
#include <Util.h>

#include <lapackpp/gmd.h>
#include <lapackpp/gmi.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>

using namespace std;

int main(int argc,char* argv[])
{
    double pe,mse;
    uint iChannelOrder,iSNR;
    int d;
    char buffer[SPRINTF_BUFFER];

    // GLOBAL PARAMETERS
    int nFrames = 2;
    int L=3,N=2,K=30;
    int trainSeqLength = 10;
    int nParticles = 10;
    double resamplingRatio = 0.9;
    char outputFileName[HOSTNAME_LENGTH+4] = "res_";
    int preambleLength = 10;

    // - ONE CHANNEL ORDER SYSTEM
    int m = 2;

    // - ONE CHANNEL ORDER PER ANTENNA SYSTEM
    vector<int> antennasChannelOrders(N);
    antennasChannelOrders[0] = 1;
    antennasChannelOrders[1] = 3;

    // SNRs to be processed
    vector<int> SNRs;
//     SNRs.push_back(3);SNRs.push_back(6);SNRs.push_back(9);SNRs.push_back(12);SNRs.push_back(15);
	SNRs.push_back(12);

    // AR process parameters
    vector<double> ARcoefficients(1);
    ARcoefficients[0] = 0.99999;
    double ARvariance=0.0001;

    // channel parameters
    double channelMean=0.0,channelVariance=1.0;

    // channel estimator parameters
	double firstSampledChannelMatrixVariance = 0.0;
    double subsequentSampledChannelMatricesVariance = 0.0;
    double forgettingFactor = 0.99;
    double muLMS = 0.05;

	// unknown channel order
	vector<int> candidateChannelOrders;
// 	candidateChannelOrders.push_back(2);candidateChannelOrders.push_back(3);candidateChannelOrders.push_back(4);

// 	candidateChannelOrders.push_back(2);candidateChannelOrders.push_back(3);candidateChannelOrders.push_back(5);candidateChannelOrders.push_back(8);candidateChannelOrders.push_back(10);

	candidateChannelOrders.push_back(2);candidateChannelOrders.push_back(3);candidateChannelOrders.push_back(4);candidateChannelOrders.push_back(5);candidateChannelOrders.push_back(6);candidateChannelOrders.push_back(7);candidateChannelOrders.push_back(8);

	// initial channel estimation for every channel order
	vector<tMatrix> initialChannelMatrixEstimations(candidateChannelOrders.size());
	for(iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
		initialChannelMatrixEstimations[iChannelOrder] = LaGenMatDouble::zeros(L,N*candidateChannelOrders[iChannelOrder]);

    // linear detectors parameters
    double forgettingFactorDetector = 0.95;

    // alphabet is defined
    vector<vector<tBit> > secuenciasBits(2,vector<tBit>(1));
    secuenciasBits[0][0] = 0; secuenciasBits[1][0] = 1;
    vector<tSymbol> simbolos(2);
    simbolos[0] = -1; simbolos[1] = 1;
    Alphabet pam2(1,2,secuenciasBits,simbolos);

	// -----------------------------------------------------------------------------

	// a vector that will contain the names of the algorithms
	vector<string> algorithmsNames;

	// host name is concatenated into the file name
	char hostname[HOSTNAME_LENGTH];
	gethostname(hostname,HOSTNAME_LENGTH);
	strcat(outputFileName,hostname);

	// get present time of the system
	time_t presentTime;
	time(&presentTime);
	char presentTimeString[DATE_LENGTH];
	ctime_r(&presentTime,presentTimeString);
	presentTimeString[strlen(presentTimeString)-1] = '\0';
	for(int i=strlen(presentTimeString)-1;i>=0;i--)
		if(presentTimeString[i]==' ')
			presentTimeString[i]='_';

	// it is concatenated into the file name
	strcat(outputFileName,"_");
	strcat(outputFileName,presentTimeString);

    // a specific preamble is generated...
    tMatrix preamble(N,preambleLength);
    preamble = -1.0;

    tRange rAllSymbolRows(0,N-1);
    tRange rTrainingSequenceSymbolVectors(preamble.cols(),preamble.cols()+trainSeqLength-1);

    // the last "d" observations are used only for smoothing and don't result in any detected vector
    tRange rSymbolVectorsToComputePe(preamble.cols()+trainSeqLength,K+preamble.cols()-1);
//     tRange rSymbolVectorsToComputePe(100,K+preamble.cols()-1);

	// vectors of channel estimators and linear detectors for unknown channel order algorithms
	vector<ChannelMatrixEstimator *> kalmanChannelEstimators;
	vector<ChannelMatrixEstimator *> RLSchannelEstimators;
	vector<LinearDetector *> RMMSElinearDetectors;
	for(iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		kalmanChannelEstimators.push_back(new KalmanEstimator(initialChannelMatrixEstimations[iChannelOrder],N,ARcoefficients[0],ARvariance));

		RLSchannelEstimators.push_back(new RLSEstimator(initialChannelMatrixEstimations[iChannelOrder],N,forgettingFactor));

		RMMSElinearDetectors.push_back(new RMMSEDetector(L*candidateChannelOrders[iChannelOrder],N*(2*candidateChannelOrders[iChannelOrder]-1),pam2.Variance(),forgettingFactorDetector,N*candidateChannelOrders[iChannelOrder]));
	}


	// channel order estimators
	APPbasedChannelOrderEstimator *channelOrderEstimator= new APPbasedChannelOrderEstimator(N,preamble,candidateChannelOrders,initialChannelMatrixEstimations,ARcoefficients[0]);

	// the maximum of the candidate channel orders is computed
	int maxCandidateOrder = candidateChannelOrders[0];
	for(iChannelOrder=1;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
		if(candidateChannelOrders[iChannelOrder]>maxCandidateOrder)
			maxCandidateOrder = candidateChannelOrders[iChannelOrder];

	// the algorithms with the higher smoothing lag require
	int nSmoothingInstants = maxCandidateOrder-1;
// 	int nSmoothingInstants = 150;

    // always the same resampling criterion and algorithms
    ResamplingCriterion criterioRemuestreo(resamplingRatio);
    StdResamplingAlgorithm algoritmoRemuestreo(criterioRemuestreo);

    // resampling algorithm for the case of unknown channel order
    ByChannelOrderResamplingAlgorithm unknownChannelOrderResamplingAlgorithm(criterioRemuestreo);

    // matrices for results
    tMatrix overallPeMatrix;
    tMatrix overallMseMatrix;

    #ifdef CHANNELORDERSAPP_SAVING
    	vector<tMatrix> channelOrdersAPPs(SNRs.size(),LaGenMatDouble::zeros(candidateChannelOrders.size(),K));
    #endif

    vector<tMatrix> overallPeTimeEvolution(SNRs.size());
    vector<LaGenMatInt> overallErrorsNumberTimeEvolution(SNRs.size());
    tMatrix channelOrderAPPsAfterTrainingSequence = LaGenMatDouble::zeros(candidateChannelOrders.size(),SNRs.size());

    // we don't want the same bits to be generated over and over
	#ifdef RANDOM_SEED
		Random bitsRandomGenerator;
	#else
    	Random bitsRandomGenerator(0);
    #endif

    for(int iFrame=0;iFrame<nFrames;iFrame++)
    {
        // bits are generated ...
        Bits bitsTransmitir(N,K+nSmoothingInstants,bitsRandomGenerator);

        // ... and then modulated by means of the alphabet
        tMatrix simbolosTransmitir = Modulator::Modulate(bitsTransmitir,pam2);

        // the preamble is set before the symbols to be transmitted
        simbolosTransmitir = Util::Append(preamble,simbolosTransmitir);

        tMatrix trainingSequence = simbolosTransmitir(rAllSymbolRows,rTrainingSequenceSymbolVectors);

        // an AR channel is generated
        ARchannel canal(N,L,m,simbolosTransmitir.cols(),channelMean,channelVariance,ARcoefficients,ARvariance);

// 		ARoneChannelOrderPerTransmitAtennaMIMOChannel canal(N,L,simbolosTransmitir.cols(),antennasChannelOrders,channelMean,channelVariance,ARcoefficients,ARvariance);

		// "m" and "d" are obtained from the just built channel object and...
		m = canal.EffectiveMemory();
		d = m-1;

		// ...according to that m, channel estimators are constructed...
		tMatrix initialChannelEstimation = LaGenMatDouble::zeros(L,N*m);

		KalmanEstimator kalmanEstimator(initialChannelEstimation,N,ARcoefficients[0],ARvariance);
	    RLSEstimator rlsEstimator(initialChannelEstimation,N,forgettingFactor);
		LMSEstimator lmsEstimator(initialChannelEstimation,N,muLMS);

		// wrapped channel estimators
// 		OneChannelOrderPerTransmitAtennaWrapperEstimator rlsWrapper(initialChannelEstimation,N,antennasChannelOrders,new RLSEstimator(OneChannelOrderPerTransmitAtennaMIMOChannel::WithZerosMatrixToWithoutZerosMatrix(initialChannelEstimation,N,antennasChannelOrders),N,forgettingFactor));

// 		OneChannelOrderPerTransmitAtennaWrapperEstimator lmsWrapper(initialChannelEstimation,N,antennasChannelOrders,new LMSEstimator(OneChannelOrderPerTransmitAtennaMIMOChannel::WithZerosMatrixToWithoutZerosMatrix(initialChannelEstimation,N,antennasChannelOrders),N,muLMS));

// 		OneChannelOrderPerTransmitAtennaWrapperEstimator kalmanWrapper(initialChannelEstimation,N,antennasChannelOrders,new KalmanEstimator(OneChannelOrderPerTransmitAtennaMIMOChannel::WithZerosMatrixToWithoutZerosMatrix(initialChannelEstimation,N,antennasChannelOrders),N,ARcoefficients[0],ARvariance));

		// ...and linear detectors
		RMMSEDetector rmmseDetector(L*(d+1),N*(m+d),pam2.Variance(),forgettingFactorDetector,N*(d+1));
// 		MMSEDetector MMSEdetector(L*(d+1),N*(m+d),pam2.Variance(),N*(d+1));

		// noise is generated according to the channel
		ChannelDependentNoise ruido(&canal);

		// absence of noise
		NullNoise ruidoNulo(canal.Nr(),canal.Length());

        for(iSNR=0;iSNR<SNRs.size();iSNR++)
        {
            cout << "SNR = " << SNRs[iSNR] << " [Trama " << iFrame << "]..." << endl;

			// noise SNR is set
            ruido.SetSNR(SNRs[iSNR],pam2.Variance());

            // transmission
            tMatrix observaciones = canal.Transmit(simbolosTransmitir,ruido);

            // algorithms are created
            vector<Algorithm *> algorithms;

            // ----------------------- ALGORITHMS TO RUN ----------------------------

//             algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt",pam2,L,N,K+preamble.cols(),m,&kalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo));

//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("LMS-D-SIS",pam2,L,N,K+preamble.cols(),m,&lmsEstimator,&rmmseDetector,preamble,d,nParticles,algoritmoRemuestreo,ARcoefficients[0],firstSampledChannelMatrixVariance,subsequentSampledChannelMatricesVariance));

//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS",pam2,L,N,K+preamble.cols(),m,&rlsEstimator,&rmmseDetector,preamble,d,nParticles,algoritmoRemuestreo,ARcoefficients[0],firstSampledChannelMatrixVariance,subsequentSampledChannelMatricesVariance));

//             algorithms.push_back(new ViterbiAlgorithm("Viterbi",pam2,L,N,K+preamble.cols(),canal,preamble,d));

//             algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimator("Kalman Filter (Known Symbols)",pam2,L,N,K+preamble.cols(),m,&kalmanEstimator,preamble,simbolosTransmitir));

							// -------- One channel order per antenna ------
//             algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt (one channel order per antenna)",pam2,L,N,K+preamble.cols(),m,&kalmanWrapper,preamble,d,nParticles,algoritmoRemuestreo));
//
//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("LMS-D-SIS (one channel order per antenna)",pam2,L,N,K+preamble.cols(),m,&lmsWrapper,&rmmseDetector,preamble,d,nParticles,algoritmoRemuestreo,ARcoefficients[0],firstSampledChannelMatrixVariance,subsequentSampledChannelMatricesVariance));
//
//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS (one channel order per antenna)",pam2,L,N,K+preamble.cols(),m,&rlsWrapper,&rmmseDetector,preamble,d,nParticles,algoritmoRemuestreo,ARcoefficients[0],firstSampledChannelMatrixVariance,subsequentSampledChannelMatricesVariance));
//
//             algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimator("Kalman Filter (Known Symbols) (one channel order per antenna)",pam2,L,N,K+preamble.cols(),m,&kalmanWrapper,preamble,simbolosTransmitir));
							// ---------------------------------------------

//             algorithms.push_back(new ISIS("ISIS",pam2,L,N,K+preamble.cols(),kalmanChannelEstimators,preamble,preamble.cols(),d,nParticles,&algoritmoRemuestreo,canal,simbolosTransmitir));

            algorithms.push_back(new USIS("USIS",pam2,L,N,K+preamble.cols(),RLSchannelEstimators,RMMSElinearDetectors,preamble,preamble.cols(),d,nParticles,&algoritmoRemuestreo,channelOrderEstimator,ARcoefficients[0],firstSampledChannelMatrixVariance,subsequentSampledChannelMatricesVariance/*,canal,simbolosTransmitir*/));

			// the RLS algorithm considering all posible channel orders
// 			for(iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
// 			{
// 				// the channel order (int) is converted to char *
// 				sprintf(buffer,"%d",candidateChannelOrders[iChannelOrder]);
//
// 				algorithms.push_back(new LinearFilterBasedSMCAlgorithm(string("Filtro lineal RLS suponiendo m = ") + string(buffer),pam2,L,N,K+preamble.cols(),candidateChannelOrders[iChannelOrder],RLSchannelEstimators[iChannelOrder],RMMSElinearDetectors[iChannelOrder],preamble,candidateChannelOrders[iChannelOrder]-1,nParticles,algoritmoRemuestreo,ARcoefficients[0],firstSampledChannelMatrixVariance,subsequentSampledChannelMatricesVariance));
// 			}

			// ---------------------------------------------------------------------------------

            // here the number of algoriths is known. So, the first iteration:
            if(iFrame==0 && iSNR==0)
            {
                overallPeMatrix.resize(SNRs.size(),algorithms.size());
                overallPeMatrix = 0.0;

                overallMseMatrix.resize(SNRs.size(),algorithms.size());
                overallMseMatrix = 0.0;

				// Pe evolution
				for(uint i=0;i<SNRs.size();i++)
				{
					overallPeTimeEvolution[i] = tMatrix(algorithms.size(),K);
					overallErrorsNumberTimeEvolution[i] = LaGenMatInt::zeros(algorithms.size(),K);
				}

				// we fill the vector with the names of the algorithms
				for(uint iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
				{
					algorithmsNames.push_back(algorithms[iAlgorithm]->GetName());
				}
            }

			#ifdef DEBUG
				vector<ChannelMatrixEstimator *> RLSchannelEstimatorsCopia = RLSchannelEstimators;
				vector<vector<tMatrix> > matricesEstimadas(candidateChannelOrders.size());
				for(iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
				{
					RLSchannelEstimatorsCopia[iChannelOrder] = RLSchannelEstimators[iChannelOrder]->Clone();
					for(int iTiempo=preambleLength;iTiempo<preambleLength+trainSeqLength;iTiempo++)
					{
						matricesEstimadas[iChannelOrder].push_back(RLSchannelEstimatorsCopia[iChannelOrder]->NextMatrix(observaciones.col(iTiempo),simbolosTransmitir(tRange(0,N-1),tRange(iTiempo-candidateChannelOrders[iChannelOrder]+1,iTiempo)),ruido.Variances()[iTiempo]));
					}
				}
				vector<double> estimatedChanelOrderAPPs = channelOrderEstimator->Clone()->ComputeProbabilities(observaciones,matricesEstimadas,ruido.Variances(),trainingSequence);
				for(iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
               		channelOrderAPPsAfterTrainingSequence(iChannelOrder,iSNR) += estimatedChanelOrderAPPs[iChannelOrder];
               	Util::Print(estimatedChanelOrderAPPs);
			#endif

            // algorithms are executed
            for(uint iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
            {
                algorithms[iAlgorithm]->Run(observaciones,ruido.Variances(),trainingSequence);

                pe = algorithms[iAlgorithm]->SER(simbolosTransmitir(rAllSymbolRows,rSymbolVectorsToComputePe));
                mse = algorithms[iAlgorithm]->MSE(canal.Range(preamble.cols(),K+preamble.cols()-1));

                cout << algorithms[iAlgorithm]->GetName() << ": Pe = " << pe << " , MSE = " << mse << endl;

                // the error probability is accumulated
                overallPeMatrix(iSNR,iAlgorithm) += pe;

                // and the MSE
                overallMseMatrix(iSNR,iAlgorithm) += mse;

                #ifdef CHANNELORDERSAPP_SAVING
                	if(!algorithms[iAlgorithm]->GetName().compare("UCO-SIS"))
                	{
                		tMatrix channelOrderAPPsAux = (dynamic_cast <USIS *>(algorithms[iAlgorithm]))->GetChannelOrderAPPsAlongTime();

                		Util::Add(channelOrderAPPsAux,channelOrdersAPPs[iSNR],channelOrdersAPPs[iSNR]);
                	}
                #endif

                // Pe evolution
                tMatrix detectedSymbols = algorithms[iAlgorithm]->GetDetectedSymbolVectors();
                tMatrix transmittedSymbols = simbolosTransmitir(tRange(0,N-1),tRange(preambleLength,preambleLength+K-1));
                #ifdef DEBUG2
					cout << "simbolos detectados" << endl << detectedSymbols << endl;
					cout << "simbolos trasmitidos" << endl << simbolosTransmitir(tRange(0,N-1),tRange(preambleLength,preambleLength+K-1)) <<endl;
                #endif
                for(int k=0;k<K;k++)
                	for(int iUser=0;iUser<N;iUser++)
                		if(detectedSymbols(iUser,k)!=transmittedSymbols(iUser,k))
                			overallErrorsNumberTimeEvolution[iSNR](iAlgorithm,k)++;

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

		for(uint iSNR=0;iSNR<SNRs.size();iSNR++)
			for(uint i=0;i<algorithmsNames.size();i++)
				for(int j=0;j<K;j++)
					overallPeTimeEvolution[iSNR](i,j) = (double) overallErrorsNumberTimeEvolution[iSNR](i,j) / (double) (N*(iFrame+1));
		Util::MatricesVectorToStream(overallPeTimeEvolution,"peTimeEvolution",f);

		tMatrix auxChannelOrderAPPsAfterTrainingSequence = channelOrderAPPsAfterTrainingSequence;
		auxChannelOrderAPPsAfterTrainingSequence *= 1.0/(double)(iFrame+1);
		Util::MatrixToStream(auxChannelOrderAPPsAfterTrainingSequence,"channelOrderAPPsAfterTrainingSequence",f);

		#ifdef CHANNELORDERSAPP_SAVING
			vector<tMatrix> auxChannelOrdersAPPs = channelOrdersAPPs;
			for(uint i=0;i<SNRs.size();i++)
				auxChannelOrdersAPPs[i] *= 1.0/(double)(iFrame+1);
			Util::MatricesVectorToStream(auxChannelOrdersAPPs,"uco_channelOrdersAPPs",f);
		#endif

        Util::ScalarToStream(iFrame+1,"nFrames",f);

		Util::StringsVectorToStream(algorithmsNames,"algorithmsNames",f);
        Util::ScalarToStream(L,"L",f);
        Util::ScalarToStream(N,"N",f);
        Util::ScalarToStream(m,"m",f);
        Util::ScalarToStream(K,"K",f);
        Util::ScalarToStream(trainSeqLength,"trainSeqLength",f);
        Util::ScalarToStream(nParticles,"nParticles",f);
        Util::ScalarToStream(resamplingRatio,"resamplingRatio",f);
        Util::ScalarToStream(d,"d",f);
		Util::ScalarsVectorToStream(candidateChannelOrders,"candidateOrders",f);
		Util::ScalarsVectorToStream(SNRs,"SNRs",f);
		Util::ScalarToStream(forgettingFactor,"forgettingFactor",f);
		Util::ScalarToStream(muLMS,"muLMS",f);
		Util::ScalarsVectorToStream(ARcoefficients,"ARcoefficients",f);
		Util::ScalarToStream(ARvariance,"ARvariance",f);
		Util::ScalarToStream(forgettingFactorDetector,"forgettingFactorDetector",f);
		Util::MatrixToStream(preamble,"preamble",f);
		Util::ScalarToStream(channelMean,"channelMean",f);
		Util::ScalarToStream(channelVariance,"channelVariance",f);
		Util::ScalarToStream(firstSampledChannelMatrixVariance,"firstSampledChannelMatrixVariance",f);
		Util::ScalarToStream(subsequentSampledChannelMatricesVariance,"subsequentSampledChannelMatricesVariance",f);
		Util::ScalarToStream(nSmoothingInstants,"nSmoothingInstants",f);
		Util::ScalarToStream(preambleLength,"preambleLength",f);

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
