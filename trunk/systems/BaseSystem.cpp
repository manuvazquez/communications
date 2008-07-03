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
#include "BaseSystem.h"

#include <defines.h>
#include <typeinfo>
#define DATE_LENGTH 100
// #define EXPORT_REAL_DATA

// #define PRINT_PARAMETERS

using namespace std;

#ifdef EXPORT_REAL_DATA
	MIMOChannel *realChannel;
	tMatrix *realSymbols;
	Noise *realNoise;
#endif

BaseSystem::BaseSystem()
{
    // GLOBAL PARAMETERS
    nFrames =1000;
    L=3,N=2,K=300;
    m = 3;
    d = m - 1;
    trainSeqLength = 10;
    sprintf(outputFileName,"res_");
    preambleLength = 10;

    SNRs.push_back(3);SNRs.push_back(6);SNRs.push_back(9);SNRs.push_back(12);SNRs.push_back(15);
//     SNRs.push_back(3);SNRs.push_back(6);
//     SNRs.push_back(9);SNRs.push_back(12);SNRs.push_back(15);
//     SNRs.push_back(3);

    // BER and MSE computing
    BERwindowStart = trainSeqLength;
    BERwindowStart = K*3/10;
    MSEwindowStart = 0;
    MSEwindowStart = K*9/10;

    // alphabet is defined
    vector<vector<tBit> > secuenciasBits(2,vector<tBit>(1));
    secuenciasBits[0][0] = 0; secuenciasBits[1][0] = 1;
    vector<tSymbol> simbolos(2);
    simbolos[0] = -1; simbolos[1] = 1;
    alphabet = new Alphabet(1,2,secuenciasBits,simbolos);

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
    preamble = tMatrix(N,preambleLength);
    preamble = -1.0;

    // the algorithms with the higher smoothing lag require
    nSmoothingSymbolsVectors = 10;
    nSmoothingBitsVectors = nSmoothingSymbolsVectors*alphabet->NbitsBySymbol();


    // ambiguity resolution
    uint *firstPermutation = new uint[N];
    for(int i=0;i<N;i++) firstPermutation[i] = i;
    permutations = Util::Permutations(firstPermutation,N);
    delete[] firstPermutation;

    peMatrices.reserve(nFrames);
    MSEMatrices.reserve(nFrames);

    overallPeTimeEvolution.resize(SNRs.size());
    overallErrorsNumberTimeEvolution.resize(SNRs.size());

    mainSeeds.reserve(nFrames);
    statUtilSeeds.reserve(nFrames);
    beforeRunStatUtilSeeds.reserve(nFrames);

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    presentFrameMSEtimeEvolution.resize(SNRs.size());
    MSEtimeEvolution.reserve(nFrames);
#endif

    // we don't want the same bits to be generated over and over
#ifndef RANDOM_SEED
        randomGenerator.setSeed(0);
#endif

	channel = NULL;
	powerProfile = NULL;
}

BaseSystem::~BaseSystem()
{
	delete alphabet;
}

void BaseSystem::Simulate()
{
  tRange rAll;
    // for repeating simulations
//     randomGenerator.setSeed();
//     StatUtil::GetRandomGenerator().setSeed();

	int iFrame = 0;
	while((iFrame<nFrames) && (!__done))
    {

        // the seeds are kept for saving later
        mainSeeds.push_back(randomGenerator.getSeed());
        statUtilSeeds.push_back(StatUtil::GetRandomGenerator().getSeed());

        // bits are generated ...
        Bits bits(N,K+nSmoothingBitsVectors,randomGenerator);

        // ... and then modulated by means of the alphabet
        symbols = Modulator::Modulate(bits,*alphabet);

        // the preamble is set before the symbols to be transmitted
        symbols = Util::Append(preamble,symbols);

        // all the above symbols must be processed except those generated due to the smoothing
        lastSymbolVectorInstant = symbols.cols() - nSmoothingSymbolsVectors;

        BuildChannel();

#ifdef PRINT_PARAMETERS
	    cout << "lastSymbolVectorInstant = " << lastSymbolVectorInstant << endl;
#endif

        // noise is generated according to the channel
// 	    ruido = new NullNoise(L,channel->Length());
        ruido = new ChannelDependentNoise(channel);
// 	    ruido = new PowerProfileDependentNoise(L,channel->Length(),*powerProfile);

#ifdef EXPORT_REAL_DATA
            realSymbols = &symbols;
            realChannel = channel;
            realNoise = ruido;
#endif

        for(iSNR=0;iSNR<SNRs.size();iSNR++)
        {
            cout << "SNR = " << SNRs[iSNR] << " [Frame " << iFrame << "]..." << endl;

            // noise SNR is set
            ruido->SetSNR(SNRs[iSNR],alphabet->Variance());

            // transmission
            observaciones = channel->Transmit(symbols,*ruido);

             AddAlgorithms();

            // here the number of algoriths is known. So, the first iteration:
            if(iFrame==0 && iSNR==0)
                OnlyOnce();

            // algorithms are executed
            for(uint iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
            {
                // in order to repeat a concrete simulation...
//                 StatUtil::GetRandomGenerator().setSeed();

                // the seed kept by the class StatUtil is saved
                presentFrameStatUtilSeeds(iSNR,iAlgorithm) = StatUtil::GetRandomGenerator().getSeed();

                algorithms[iAlgorithm]->Run(observaciones,ruido->Variances(),symbols(rAll,tRange(preambleLength,preambleLength+trainSeqLength-1)));
//                 algorithms[iAlgorithm]->Run(observaciones,ruido->Variances());

                detectedSymbols = algorithms[iAlgorithm]->GetDetectedSymbolVectors();

                pe = TransmissionUtil::ComputeBERsolvingAmbiguity(bits,BERwindowStart,K,Demodulator::Demodulate(detectedSymbols,*alphabet),BERwindowStart,K,permutations);

                BeforeEndingAlgorithm(iAlgorithm);

                delete algorithms[iAlgorithm];
            }

            algorithms.clear();
        } // for(int iSNR=0;iSNR<SNRs.size();iSNR++)

        f.open(outputFileName,ofstream::trunc);
        BeforeEndingFrame(iFrame);
        f.close();

        // ---------------------------------------------------------

	    iFrame++;

	    delete channel;
	    delete ruido;
    } // while((iFrame<nFrames) && (!done))

	overallPeMatrix *= 1.0/iFrame;
	overallMseMatrix *= 1.0/iFrame;

    cout << "Overall SER:" << endl;
    Util::Print(overallPeMatrix);

    cout << "Overall MSE:" << endl;
    Util::Print(overallMseMatrix);
}

void BaseSystem::OnlyOnce()
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

    // we fill the vector with the names of the algorithms
    for(uint iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
        algorithmsNames.push_back(algorithms[iAlgorithm]->GetName());

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    for(uint i=0;i<SNRs.size();i++)
        presentFrameMSEtimeEvolution[i] = tMatrix(algorithms.size(),K);
#endif

    presentFrameStatUtilSeeds = LaGenMatLongInt(SNRs.size(),algorithms.size());
}

void BaseSystem::BeforeEndingFrame(int iFrame)
{
    // pe
    peMatrices.push_back(presentFramePe);
    Util::MatricesVectorToOctaveFileStream(peMatrices,"pe",f);

    // MSE
    MSEMatrices.push_back(presentFrameMSE);
    Util::MatricesVectorToOctaveFileStream(MSEMatrices,"mse",f);

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    MSEtimeEvolution.push_back(presentFrameMSEtimeEvolution);
    Util::MatricesVectoresVectorToOctaveFileStream(MSEtimeEvolution,"MSEtimeEvolution",f);
#endif

    // seeds just before the run of the algorithms
    beforeRunStatUtilSeeds.push_back(presentFrameStatUtilSeeds);
    Util::MatricesVectorToOctaveFileStream(beforeRunStatUtilSeeds,"beforeRunStatUtilSeeds",f);

//     for(uint iSNR=0;iSNR<SNRs.size();iSNR++)
//         for(uint i=0;i<algorithmsNames.size();i++)
//             for(int j=0;j<K;j++)
//                 overallPeTimeEvolution[iSNR](i,j) = (double) overallErrorsNumberTimeEvolution[iSNR](i,j) / (double) (N*(iFrame+1));
//     Util::MatricesVectorToOctaveFileStream(overallPeTimeEvolution,"peTimeEvolution",f);

    Util::ScalarToOctaveFileStream(iFrame+1,"nFrames",f);

    Util::StringsVectorToOctaveFileStream(algorithmsNames,"algorithmsNames",f);
    Util::ScalarToOctaveFileStream(L,"L",f);
    Util::ScalarToOctaveFileStream(N,"N",f);
    Util::ScalarToOctaveFileStream(m,"m",f);
    Util::ScalarToOctaveFileStream(K,"K",f);
    Util::ScalarToOctaveFileStream(trainSeqLength,"trainSeqLength",f);
    Util::ScalarToOctaveFileStream(d,"d",f);
    Util::ScalarToOctaveFileStream(BERwindowStart,"BERwindowStart",f);
    Util::ScalarToOctaveFileStream(MSEwindowStart,"MSEwindowStart",f);
    Util::ScalarsVectorToOctaveFileStream(SNRs,"SNRs",f);
    Util::MatrixToOctaveFileStream(preamble,"preamble",f);
    Util::ScalarToOctaveFileStream(nSmoothingBitsVectors,"nSmoothingBitsVectors",f);
    Util::ScalarToOctaveFileStream(preambleLength,"preambleLength",f);
    Util::ScalarsVectorToOctaveFileStream(mainSeeds,"mainSeeds",f);
    Util::ScalarsVectorToOctaveFileStream(statUtilSeeds,"statUtilSeeds",f);
//     Util::MatricesVectorToOctaveFileStream(channel->Range(preambleLength,lastSymbolVectorInstant),"channel",f);
	Util::StringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*channel).name())),"channelClass",f);
	Util::StringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*ruido).name())),"ruidoClass",f);
    Util::StringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*this).name())),"systemClass",f);

	if(powerProfile!=NULL)
	{
		Util::ScalarsVectorToOctaveFileStream(powerProfile->TapsAmplitudes(),"powerProfileVariances",f);
		Util::StringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*powerProfile).name())),"powerProfileClass",f);
	}
}

void BaseSystem::BeforeEndingAlgorithm(int iAlgorithm)
{
    mse = algorithms[iAlgorithm]->MSE(channel->Range(preambleLength+MSEwindowStart,lastSymbolVectorInstant-1));

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    tVector mseAlongTime = TransmissionUtil::MSEalongTime(algorithms[iAlgorithm]->GetEstimatedChannelMatrices(),0,K-1,channel->Range(preambleLength,preambleLength+K-1),0,K-1);
    for(int ik=0;ik<K;ik++)
        presentFrameMSEtimeEvolution[iSNR](iAlgorithm,ik) = mseAlongTime(ik);
#endif

    cout << algorithms[iAlgorithm]->GetName() << ": Pe = " << pe << " , MSE = " << mse << endl;

    // the error probability is accumulated
    overallPeMatrix(iSNR,iAlgorithm) += pe;
    presentFramePe(iSNR,iAlgorithm) = pe;

    // and the MSE
    overallMseMatrix(iSNR,iAlgorithm) += mse;
    presentFrameMSE(iSNR,iAlgorithm) = mse;

    // Pe evolution
    tMatrix transmittedSymbols = symbols(tRange(0,N-1),tRange(preambleLength,preambleLength+K-1));

    if(detectedSymbols.rows()!=0)
    {
        for(int k=0;k<K;k++)
            for(int iUser=0;iUser<N;iUser++)
                if(detectedSymbols(iUser,k)!=transmittedSymbols(iUser,k))
                    overallErrorsNumberTimeEvolution[iSNR](iAlgorithm,k)++;
    }
}
