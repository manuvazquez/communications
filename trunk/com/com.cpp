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
#define SPRINTF_BUFFER 50

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <unistd.h>
#include <string>
#include <math.h>
#include <sys/time.h>
#include <algorithm>

#include <types.h>
#include <defines.h>
#include <Alphabet.h>
#include <Bits.h>
#include <ARprocess.h>

#include <StatUtil.h>
#include <Util.h>

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
#include <LinearFilterBasedMKFAlgorithm.h>
#include <ViterbiAlgorithm.h>
#include <KnownSymbolsKalmanBasedChannelEstimator.h>
#include <UnknownChannelOrderAlgorithm.h>
#include <ISIS.h>
#include <USIS.h>
#include <USIS2SISAlgorithm.h>
#include <MaximumProbabilityCriterion.h>
#include <UniformRelatedCriterion.h>
#include <PSPAlgorithm.h>
#include <StochasticPSPAlgorithm.h>
#include <PSPBasedSMCAlgorithm.h>

#include <Particle.h>
#include <ParticleWithChannelEstimation.h>
#include <ResamplingCriterion.h>
#include <MultinomialResamplingAlgorithm.h>
#include <ResidualResamplingAlgorithm.h>
#include <WithThresholdResamplingAlgorithmWrapper.h>
#include <WithoutReplacementResamplingAlgorithm.h>
#include <BestParticlesResamplingAlgorithm.h>

#include <lapackpp/gmd.h>
#include <lapackpp/gmi.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>
#include <lapackpp/sybmd.h>
#include <lapackpp/sybfd.h>

using namespace std;

double ComputeBER(const Bits &bits1,int from1,int to1,const Bits &bits2,int from2,int to2);
double ComputeBERsolvingAmbiguity(const Bits &sourceBits,int from1,int to1,const Bits &detectedBits,int from2,int to2,vector<vector<uint> > permutations);
void BERComputingChecks(const Bits &sourceBits,int from1,int to1,const Bits &detectedBits,int from2,int to2);

#ifdef EXPORT_REAL_DATA
	MIMOChannel *realChannel;
	tMatrix *realSymbols;
	Noise *realNoise;
#endif

int main(int argc,char* argv[])
{
// 	tMatrix matriz(2,2);
// 	matriz(0,0) = 1;matriz(1,0) = 2;matriz(0,1) = 3;matriz(1,1) = 4;
//
// 	vector<tMatrix> vectorMatrices(4,matriz);
// 	for(uint i=0;i<vectorMatrices.size();i++)
// 		vectorMatrices[i] += double(i*10);
//
// 	vector<vector<tMatrix> > vectorVectoresMatrices(3,vectorMatrices);
// 	for(uint i=0;i<vectorVectoresMatrices.size();i++)
// 		for(uint j=0;j<vectorVectoresMatrices[i].size();j++)
// 			vectorVectoresMatrices[i][j] += double(i*100);
//
// 	vector<vector<vector<tMatrix> > > vectorVectoresVectoresMatrices(2,vectorVectoresMatrices);
// 	for(uint i=0;i<vectorVectoresVectoresMatrices.size();i++)
// 		for(uint j=0;j<vectorVectoresVectoresMatrices[i].size();j++)
// 			for(uint k=0;k<vectorVectoresVectoresMatrices[i][j].size();k++)
// 				vectorVectoresVectoresMatrices[i][j][k] += double(i*1000);
//
// 	ofstream f2("venga",ofstream::trunc);
// 	Util::MatricesVectoresVectoresVectorToStream(vectorVectoresVectoresMatrices,"si",f2);
// 	f2.close();
// 	exit(0);


    double pe,mse;
    uint iChannelOrder,iSNR;
    int d,lastSymbolVectorInstant;

    // GLOBAL PARAMETERS
    int nFrames = 2;
    int L=3,N=2,K=300;
    int trainSeqLength = 30;
    int nParticles = 30;
    double resamplingRatio = 0.9;
    char outputFileName[HOSTNAME_LENGTH+4] = "res_";
    int preambleLength = 10;

    // BER and MSE computing
    int BERwindowStart = trainSeqLength;
    int MSEwindowStart = 0;

    // PSP
    int nSurvivors = 2;
    bool adjustParticlesNumberFromSurvivors = false;

    // - ONE CHANNEL ORDER SYSTEM
    int m = 4;

    // - ONE CHANNEL ORDER PER ANTENNA SYSTEM
    vector<int> antennasChannelOrders(N);
    antennasChannelOrders[0] = 1;
    antennasChannelOrders[1] = 3;

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
	double firstSampledChannelMatrixVariance = 0.0;
    double forgettingFactor = 0.99;
    double muLMS = 0.05;

	// unknown channel order
	vector<int> candidateChannelOrders;
	candidateChannelOrders.push_back(2);candidateChannelOrders.push_back(3);candidateChannelOrders.push_back(4);candidateChannelOrders.push_back(5);candidateChannelOrders.push_back(6);

	if(find(candidateChannelOrders.begin(),candidateChannelOrders.end(),m)==candidateChannelOrders.end())
		throw RuntimeException("The memory of the channel is not one of the possible candidates.");

	// initial channel estimation for every channel order
	vector<tMatrix> channelOrderCoefficientsMeans(candidateChannelOrders.size());
	vector<tMatrix> channelOrderCoefficientsVariances(candidateChannelOrders.size());
	for(iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		channelOrderCoefficientsMeans[iChannelOrder] = LaGenMatDouble::ones(L,N*candidateChannelOrders[iChannelOrder]);
		channelOrderCoefficientsMeans[iChannelOrder] *= channelMean;

		channelOrderCoefficientsVariances[iChannelOrder] = LaGenMatDouble::ones(L,N*candidateChannelOrders[iChannelOrder]);
		channelOrderCoefficientsVariances[iChannelOrder] *= channelVariance;
	}

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

    // some useful ranges
    tRange rAllSymbolRows(0,N-1);
    tRange rTrainingSequenceSymbolVectors(preambleLength,preambleLength+trainSeqLength-1);
    tRange rAllObservationsRows(0,L-1),rLastNchannelMatrixColumns(N*m-N,N*m-1);

    // variances for generating the channel coefficients
    vector<double> subChannelMatrixVariances(m);
    tMatrix channelCoefficientsVariances = LaGenMatDouble::ones(L,N*m);
    for(int i=0;i<m;i++)
    {
        subChannelMatrixVariances[i] = channelVariance;
//         subChannelMatrixVariances[i] = exp((i+1-m));
        channelCoefficientsVariances(rAllObservationsRows,tRange(i*N,i*N+N-1)) *= subChannelMatrixVariances[i];
    }

	// vectors of channel estimators and linear detectors for unknown channel order algorithms
	vector<ChannelMatrixEstimator *> kalmanChannelEstimators;
	vector<ChannelMatrixEstimator *> RLSchannelEstimators;
	vector<LinearDetector *> RMMSElinearDetectors;
	for(iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		kalmanChannelEstimators.push_back(new KalmanEstimator(channelOrderCoefficientsMeans[iChannelOrder],N,ARcoefficients[0],ARvariance));

		RLSchannelEstimators.push_back(new RLSEstimator(channelOrderCoefficientsMeans[iChannelOrder],N,forgettingFactor));

		RMMSElinearDetectors.push_back(new RMMSEDetector(L*candidateChannelOrders[iChannelOrder],N*(2*candidateChannelOrders[iChannelOrder]-1),pam2.Variance(),forgettingFactorDetector,N*candidateChannelOrders[iChannelOrder]));
	}

	// channel order estimators
	APPbasedChannelOrderEstimator *channelOrderEstimator= new APPbasedChannelOrderEstimator(N,preamble,candidateChannelOrders,channelOrderCoefficientsMeans,ARcoefficients[0]);

	// the maximum of the candidate channel orders is computed
	int maxCandidateOrder = candidateChannelOrders[0];
	for(iChannelOrder=1;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
		if(candidateChannelOrders[iChannelOrder]>maxCandidateOrder)
			maxCandidateOrder = candidateChannelOrders[iChannelOrder];

	// the algorithms with the higher smoothing lag require
	int nSmoothingSymbolsVectors = maxCandidateOrder-1;

	int nSmoothingBitsVectors = nSmoothingSymbolsVectors*pam2.NbitsBySymbol();

	// PSP
	if(adjustParticlesNumberFromSurvivors)
	{
		nParticles = (int)pow((double)pam2.Length(),N*(m-1))*nSurvivors;
        cout << "Number of particles adjusted to " << nParticles << endl;
    }

    // always the same resampling criterion and algorithms
    ResamplingCriterion criterioRemuestreo(resamplingRatio);
    MultinomialResamplingAlgorithm algoritmoRemuestreo(criterioRemuestreo);
    ResidualResamplingAlgorithm residualResampling(criterioRemuestreo);
    WithThresholdResamplingAlgorithmWrapper residualResamplingWithThreshold(new ResidualResamplingAlgorithm(criterioRemuestreo),0.2);
    WithoutReplacementResamplingAlgorithm withoutReplacementResampling(criterioRemuestreo);
    BestParticlesResamplingAlgorithm bestParticlesResampling(criterioRemuestreo);


	vector<double> resamplingRates;
	resamplingRates.push_back(0.001);resamplingRates.push_back(0.05);resamplingRates.push_back(0.1);
	resamplingRates.push_back(0.3);resamplingRates.push_back(0.5);resamplingRates.push_back(0.9);

	vector<ResidualResamplingAlgorithm *> resamplingAlgorithms;
	for(uint iResamplingAlgorithm=0;iResamplingAlgorithm<resamplingRates.size();iResamplingAlgorithm++)
		resamplingAlgorithms.push_back(new ResidualResamplingAlgorithm(ResamplingCriterion(resamplingRates[iResamplingAlgorithm])));

	// USIS2SIS transition criterion(s)
    MaximumProbabilityCriterion USISmaximumProbabilityCriterion(0.8);
    UniformRelatedCriterion USISuniformRelatedCriterion(2.0);

    // ambiguity resolution
    uint *firstPermutation = new uint[N];
    for(int i=0;i<N;i++) firstPermutation[i] = i;
    vector<vector<uint> > permutations = Util::Permutations(firstPermutation,N);
    delete[] firstPermutation;

    // matrices for results
    vector<tMatrix> peMatrices, MSEMatrices;
    peMatrices.reserve(nFrames);
    MSEMatrices.reserve(nFrames);

    tMatrix overallPeMatrix,overallMseMatrix,presentFramePe,presentFrameMSE;

   	#ifdef CHANNELORDERSAPP_SAVING
    	vector<vector<vector<tMatrix> > > channelOrderAPPestimations;
    	channelOrderAPPestimations.reserve(nFrames);

    	vector<vector<tMatrix> > presentFrameChannelOrderAPPevolution;
    #endif

    vector<tMatrix> overallPeTimeEvolution(SNRs.size());
    vector<LaGenMatInt> overallErrorsNumberTimeEvolution(SNRs.size());

    // we don't want the same bits to be generated over and over
	#ifdef RANDOM_SEED
		Random bitsRandomGenerator;
	#else
    	Random bitsRandomGenerator(0);
    #endif

    for(int iFrame=0;iFrame<nFrames;iFrame++)
    {
        // bits are generated ...
        Bits bits(N,K+nSmoothingBitsVectors,bitsRandomGenerator);

        // ... and then modulated by means of the alphabet
        tMatrix symbols = Modulator::Modulate(bits,pam2);

        // the preamble is set before the symbols to be transmitted
        symbols = Util::Append(preamble,symbols);

		#ifdef EXPORT_REAL_DATA
			realSymbols = &symbols;
		#endif

		// all the above symbols must be processed except those generated due to the smoothing
        lastSymbolVectorInstant = symbols.cols() - nSmoothingSymbolsVectors;

        tMatrix trainingSequence = symbols(rAllSymbolRows,rTrainingSequenceSymbolVectors);

        // ARchannel matrix intialization
        tMatrix ARchannelInitializationMatrix(L,N*m);

        for(int i=0;i<m;i++)
            ARchannelInitializationMatrix(rAllObservationsRows,tRange(i*N,i*N+N-1)).inject(StatUtil::RandnMatrix(L,N,channelMean,subChannelMatrixVariances[i]));

        // an AR channel is generated
        ARchannel canal(N,L,m,symbols.cols(),ARchannelInitializationMatrix,ARcoefficients,ARvariance);

// 		ARoneChannelOrderPerTransmitAtennaMIMOChannel canal(N,L,symbols.cols(),antennasChannelOrders,channelMean,channelVariance,ARcoefficients,ARvariance);

		#ifdef EXPORT_REAL_DATA
			realChannel = &canal;
		#endif

		// "m" and "d" are obtained from the just built channel object ...
		m = canal.EffectiveMemory();
		d = m-1;

		// ... and according to that m, channel estimators are constructed...
		tMatrix initialChannelEstimation = LaGenMatDouble::zeros(L,N*m);

// 		KalmanEstimator kalmanEstimator(initialChannelEstimation,N,ARcoefficients[0],ARvariance);
        KalmanEstimator kalmanEstimator(initialChannelEstimation,channelCoefficientsVariances,N,ARcoefficients[0],ARvariance);
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

		#ifdef EXPORT_REAL_DATA
			realNoise = &ruido;
		#endif

        for(iSNR=0;iSNR<SNRs.size();iSNR++)
        {
            cout << "SNR = " << SNRs[iSNR] << " [Trama " << iFrame << "]..." << endl;

			// noise SNR is set
            ruido.SetSNR(SNRs[iSNR],pam2.Variance());

            // transmission
            tMatrix observaciones = canal.Transmit(symbols,ruido);

            // algorithms are created
            vector<Algorithm *> algorithms;

            // ----------------------- ALGORITHMS TO RUN ----------------------------

//             algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt",pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,preamble,d,nParticles,&algoritmoRemuestreo,initialChannelEstimation,channelCoefficientsVariances));

//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("LMS-D-SIS",pam2,L,N,lastSymbolVectorInstant,m,&lmsEstimator,&rmmseDetector,preamble,d,nParticles,&algoritmoRemuestreo,initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

            algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS",pam2,L,N,lastSymbolVectorInstant,m,&rlsEstimator,&rmmseDetector,preamble,d,nParticles,&algoritmoRemuestreo,initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS with residual resampling",pam2,L,N,lastSymbolVectorInstant,m,&rlsEstimator,&rmmseDetector,preamble,d,nParticles,&residualResampling,initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,&canal,&symbols));

//             algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF",pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,&rmmseDetector,preamble,d,nParticles,&algoritmoRemuestreo,initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

//             algorithms.push_back(new ViterbiAlgorithm("Viterbi",pam2,L,N,lastSymbolVectorInstant,canal,preamble,d));

//             algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimator("Kalman Filter (Known Symbols)",pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,preamble,symbols));

//             algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,preamble,d,lastSymbolVectorInstant+d,ARcoefficients[0],nSurvivors));

// 			algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm",pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,preamble,d,nParticles,&withoutReplacementResampling,initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0]));

// 			algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (best particles resampling)",pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,preamble,d,nParticles,&bestParticlesResampling,initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0]));

							// -------- One channel order per antenna ------
//             algorithms.push_back(new DSISoptAlgorithm ("D-SIS opt (one channel order per antenna)",pam2,L,N,lastSymbolVectorInstant,m,&kalmanWrapper,preamble,d,nParticles,&algoritmoRemuestreo));
//
//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("LMS-D-SIS (one channel order per antenna)",pam2,L,N,lastSymbolVectorInstant,m,&lmsWrapper,&rmmseDetector,preamble,d,nParticles,&algoritmoRemuestreo,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
//
//             algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS (one channel order per antenna)",pam2,L,N,lastSymbolVectorInstant,m,&rlsWrapper,&rmmseDetector,preamble,d,nParticles,&algoritmoRemuestreo,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
//
//             algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimator("Kalman Filter (Known Symbols) (one channel order per antenna)",pam2,L,N,lastSymbolVectorInstant,m,&kalmanWrapper,preamble,symbols));
							// ---------------------------------------------

//             algorithms.push_back(new ISIS("ISIS",pam2,L,N,lastSymbolVectorInstant,kalmanChannelEstimators,preamble,preamble.cols(),d,nParticles,&algoritmoRemuestreo,canal,symbols));

//             algorithms.push_back(new USIS("USIS",pam2,L,N,lastSymbolVectorInstant,RLSchannelEstimators,RMMSElinearDetectors,preamble,preamble.cols(),d,nParticles,&algoritmoRemuestreo,channelOrderEstimator,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance/*,canal,symbols*/));

            algorithms.push_back(new USIS2SISAlgorithm("USIS2SISAlgorithm (maximum probability criterion)",pam2,L,N,lastSymbolVectorInstant,RLSchannelEstimators,RMMSElinearDetectors,preamble,preamble.cols(),d,nParticles,&algoritmoRemuestreo,channelOrderEstimator,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,&USISmaximumProbabilityCriterion/*,canal,symbols*/));


//             algorithms.push_back(new USIS2SISAlgorithm("USIS2SISAlgorithm (uniform criterion)",pam2,L,N,lastSymbolVectorInstant,RLSchannelEstimators,RMMSElinearDetectors,preamble,preamble.cols(),d,nParticles,&algoritmoRemuestreo,channelOrderEstimator,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,&USISuniformRelatedCriterion/*,canal,symbols*/));

			// the RLS algorithm considering all posible channel orders
// 			for(iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
// 			{
//
//                 char buffer[SPRINTF_BUFFER];
//
//                 if(candidateChannelOrders[iChannelOrder]<m)
//                 	sprintf(buffer," underestimated (m=%d)",candidateChannelOrders[iChannelOrder]);
//                 else if(candidateChannelOrders[iChannelOrder]>m)
//                 	sprintf(buffer," overestimated (m=%d)",candidateChannelOrders[iChannelOrder]);
//                 else buffer[0] = '\0';
//
// 				algorithms.push_back(new LinearFilterBasedSMCAlgorithm(string("RLS-D-SIS") + string(buffer),pam2,L,N,lastSymbolVectorInstant,candidateChannelOrders[iChannelOrder],RLSchannelEstimators[iChannelOrder],RMMSElinearDetectors[iChannelOrder],preamble,candidateChannelOrders[iChannelOrder]-1,nParticles,&algoritmoRemuestreo,channelOrderCoefficientsMeans[iChannelOrder],channelOrderCoefficientsVariances[iChannelOrder],ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
// 			}

// 			the RLS algorithm considering differente resampling rates
// 			for(uint iResamplingAlgorithm=0;iResamplingAlgorithm<resamplingRates.size();iResamplingAlgorithm++)
// 			{
//                 char buffer[SPRINTF_BUFFER];
//
//                 sprintf(buffer," resampling rate = %f",resamplingRates[iResamplingAlgorithm]);
//
// 				algorithms.push_back(new LinearFilterBasedSMCAlgorithm(string("RLS-D-SIS") + string(buffer),pam2,L,N,lastSymbolVectorInstant,m,&rlsEstimator,&rmmseDetector,preamble,d,nParticles,resamplingAlgorithms[iResamplingAlgorithm],initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,&canal,&symbols));
// 			}

			// ---------------------------------------------------------------------------------

            // here the number of algoriths is known. So, the first iteration:
            if(iFrame==0 && iSNR==0)
            {
                overallPeMatrix = LaGenMatDouble::zeros(SNRs.size(),algorithms.size());
                presentFramePe = LaGenMatDouble::zeros(SNRs.size(),algorithms.size());

                overallMseMatrix = LaGenMatDouble::zeros(SNRs.size(),algorithms.size());
                presentFrameMSE = LaGenMatDouble::zeros(SNRs.size(),algorithms.size());

				// Pe evolution
				for(uint i=0;i<SNRs.size();i++)
				{
					overallPeTimeEvolution[i] = tMatrix(algorithms.size(),K);
					overallErrorsNumberTimeEvolution[i] = LaGenMatInt::zeros(algorithms.size(),K);
				}

				// the number of algorithms that perform channel order APP estimation along time
				int nAlgorithmsPerformingChannelOrderAPPestimation = 0;

				// we fill the vector with the names of the algorithms
				for(uint iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
				{
					algorithmsNames.push_back(algorithms[iAlgorithm]->GetName());

					if(algorithms[iAlgorithm]->PerformsChannelOrderAPPEstimation())
						nAlgorithmsPerformingChannelOrderAPPestimation++;
				}

				#ifdef CHANNELORDERSAPP_SAVING
					// channel order APP evolution
					presentFrameChannelOrderAPPevolution = vector<vector<tMatrix> >(nAlgorithmsPerformingChannelOrderAPPestimation,vector<tMatrix>(SNRs.size(),LaGenMatDouble::zeros(candidateChannelOrders.size(),K)));
				#endif
            }

			#ifdef CHANNELORDERSAPP_SAVING
				int iAlgorithmPerformingChannelOrderAPPestimation = 0;
			#endif

            // algorithms are executed
            for(uint iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
            {
                algorithms[iAlgorithm]->Run(observaciones,ruido.Variances(),trainingSequence);
//                 algorithms[iAlgorithm]->Run(observaciones,ruido.Variances());

                tMatrix detectedSymbols = algorithms[iAlgorithm]->GetDetectedSymbolVectors();
                vector<tMatrix> estimatedChannelMatrices = algorithms[iAlgorithm]->GetEstimatedChannelMatrices();

				pe = ComputeBERsolvingAmbiguity(bits,BERwindowStart,K,Demodulator::Demodulate(detectedSymbols,pam2),BERwindowStart,K,permutations);

                mse = algorithms[iAlgorithm]->MSE(canal.Range(preambleLength+MSEwindowStart,lastSymbolVectorInstant-1));

                cout << algorithms[iAlgorithm]->GetName() << ": Pe = " << pe << " , MSE = " << mse << endl;

                // the error probability is accumulated
                overallPeMatrix(iSNR,iAlgorithm) += pe;
                presentFramePe(iSNR,iAlgorithm) = pe;

                // and the MSE
                overallMseMatrix(iSNR,iAlgorithm) += mse;
                presentFrameMSE(iSNR,iAlgorithm) = mse;

                #ifdef CHANNELORDERSAPP_SAVING
                	if(algorithms[iAlgorithm]->PerformsChannelOrderAPPEstimation())
                	{
                		presentFrameChannelOrderAPPevolution[iAlgorithmPerformingChannelOrderAPPestimation][iSNR] = (dynamic_cast <USIS *>(algorithms[iAlgorithm]))->GetChannelOrderAPPsAlongTime();
                		iAlgorithmPerformingChannelOrderAPPestimation++;

//                 		cout << "La matriz que devuelve " << endl << (dynamic_cast <USIS *>(algorithms[iAlgorithm]))->GetChannelOrderAPPsAlongTime();
//                 		cout << "Una tecla..."; getchar();
                	}
                #endif

                // Pe evolution
                tMatrix transmittedSymbols = symbols(tRange(0,N-1),tRange(preambleLength,preambleLength+K-1));

                for(int k=0;k<K;k++)
                	for(int iUser=0;iUser<N;iUser++)
                		if(detectedSymbols(iUser,k)!=transmittedSymbols(iUser,k))
                			overallErrorsNumberTimeEvolution[iSNR](iAlgorithm,k)++;

                delete algorithms[iAlgorithm];
            }
        } // for(int iSNR=0;iSNR<SNRs.size();iSNR++)


		// ----------------- VARIABLES SAVING ----------------------
		ofstream f(outputFileName,ofstream::trunc);

		peMatrices.push_back(presentFramePe);
		Util::MatricesVectorToStream(peMatrices,"pe",f);

		MSEMatrices.push_back(presentFrameMSE);
		Util::MatricesVectorToStream(MSEMatrices,"mse",f);

		#ifdef CHANNELORDERSAPP_SAVING
			channelOrderAPPestimations.push_back(presentFrameChannelOrderAPPevolution);
			Util::MatricesVectoresVectoresVectorToStream(channelOrderAPPestimations,"channelOrderAPPevolution",f);
		#endif

		for(uint iSNR=0;iSNR<SNRs.size();iSNR++)
			for(uint i=0;i<algorithmsNames.size();i++)
				for(int j=0;j<K;j++)
					overallPeTimeEvolution[iSNR](i,j) = (double) overallErrorsNumberTimeEvolution[iSNR](i,j) / (double) (N*(iFrame+1));
		Util::MatricesVectorToStream(overallPeTimeEvolution,"peTimeEvolution",f);

        Util::ScalarToStream(iFrame+1,"nFrames",f);

		Util::StringsVectorToStream(algorithmsNames,"algorithmsNames",f);
        Util::ScalarToStream(L,"L",f);
        Util::ScalarToStream(N,"N",f);
        Util::ScalarToStream(m,"m",f);
        Util::ScalarToStream(K,"K",f);
        Util::ScalarToStream(trainSeqLength,"trainSeqLength",f);
        Util::ScalarToStream(nParticles,"nParticles",f);
        Util::ScalarToStream(nSurvivors,"nSurvivors",f);
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
		Util::ScalarToStream(nSmoothingBitsVectors,"nSmoothingBitsVectors",f);
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

    // memory is released
	for(iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
	{
		delete kalmanChannelEstimators[iChannelOrder];
		delete RLSchannelEstimators[iChannelOrder];
		delete RMMSElinearDetectors[iChannelOrder];
	}
	delete channelOrderEstimator;

	for(uint iResamplingAlgorithm=0;iResamplingAlgorithm<resamplingRates.size();iResamplingAlgorithm++)
		delete resamplingAlgorithms[iResamplingAlgorithm];
}

void BERComputingChecks(const Bits &bits1,int from1,int to1,const Bits &bits2,int from2,int to2)
{
    if((to1-from1)!=(to2-from2))
   	{
   		cout << "Range 1: " << (to1-from1) << " | " << "Range 2: " << to2-from2 << endl;
        throw RuntimeException("BERComputingChecks: comparisons range length are different.");
    }

    if(to1<from1)
        throw RuntimeException("BERComputingChecks: comparisons range are negatives.");

    if(to1>bits1.NbitsByStream() || to2>bits2.NbitsByStream() || from1<0 || from2<0)
    {
    	cout << "bits1.NbitsByStream(): " << bits1.NbitsByStream() << ",bits2.NbitsByStream(): " << bits2.NbitsByStream() << endl;
    	cout << "from1: " << from1 << ", to1: " << to1 << " | " << "from2: " << from2 << ", to2: " << to2 << endl;
        throw RuntimeException("BERComputingChecks: one or several comparison limits are wrong.");
    }

    if(bits1.Nstreams()!=bits2.Nstreams())
        throw RuntimeException("BERComputingChecks: bits objects have different number of streams.");
}

double ComputeBER(const Bits &bits1,int from1,int to1,const Bits &bits2,int from2,int to2)
{
	BERComputingChecks(bits1,from1,to1,bits2,from2,to2);

	int length = to1-from1;
    int errors = 0;

    for(int iBits1=from1,iBits2=from2;iBits1<to1;iBits1++,iBits2++)
        for(int iStream=0;iStream<bits1.Nstreams();iStream++)
            if(bits1(iStream,iBits1)!=bits2(iStream,iBits2))
                errors++;

    return (double)errors/(double)(length*bits1.Nstreams());
}

double ComputeBERsolvingAmbiguity(const Bits &sourceBits,int from1,int to1,const Bits &detectedBits,int from2,int to2,vector<vector<uint> > permutations)
{
	BERComputingChecks(sourceBits,from1,to1,detectedBits,from2,to2);

	int length = to1-from1;

	// max number of errors is length*sourceBits.Nstreams()
	int minErrors = length*sourceBits.Nstreams()+1;

    for(uint iPermut=0;iPermut<permutations.size();iPermut++)
    {
    	int errorsPermutation = 0;

        for(uint iStream=0;iStream<permutations[iPermut].size();iStream++)
        {
        	int errorsInverting=0,errorsWithoutInverting=0;

        	// without inverting bits
        	for(int iSourceStream=from1,iDetectedStream=from2;iSourceStream<to1;iSourceStream++,iDetectedStream++)
        		errorsWithoutInverting += (sourceBits(iStream,iSourceStream) != detectedBits(permutations[iPermut][iStream],iDetectedStream));

        	// inverting bits
        	for(int iSourceStream=from1,iDetectedStream=from2;iSourceStream<to1;iSourceStream++,iDetectedStream++)
        		errorsInverting += (sourceBits(iStream,iSourceStream) == detectedBits(permutations[iPermut][iStream],iDetectedStream));

        	if(errorsWithoutInverting<errorsInverting)
        		errorsPermutation += errorsWithoutInverting;
        	else
        		errorsPermutation += errorsInverting;
        }

        if(errorsPermutation<minErrors)
        	minErrors = errorsPermutation;
    }

    return (double)minErrors/(double)(length*sourceBits.Nstreams());
}
