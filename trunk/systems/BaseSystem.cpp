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
#include <string.h>

#include <SingleUserPowerProfileDependentNoise.h>

#define DATE_LENGTH 100

#define EXPORT_REAL_DATA
#define PRINT_PARAMETERS
// #define PRINT_SYMBOLS_ACCOUNTED_FOR_DETECTION

// #define DEBUG

using namespace std;

#ifdef EXPORT_REAL_DATA
    MIMOChannel *realChannel;
    MatrixXd *realSymbols;
    Noise *realNoise;
#endif

BaseSystem::BaseSystem()
{
    // GLOBAL PARAMETERS
    nFrames = 10;
    L=3,N=2,frameLength=300;
    m = 3;
    d = m - 1;
    trainSeqLength = 10;
    preambleLength = 10;
  
    // the algorithms with the higher smoothing lag require
    nSmoothingSymbolsVectors = 10;
    
// 	nFrames = 2000;
// // 	nFrames = 10;
// // 	nFrames = 1;
// // 	nFrames = 200;
// //     L=3,N=2,frameLength=300;
//     L=7,N=3,frameLength=10;	
// //     L=7,N=1,frameLength=10;
// //     L=7,N=3,frameLength=300;	
//     m = 1;
//     d = m - 1;
//     trainSeqLength = 0;
//     preambleLength = 0;
//     
//     // the algorithms with the higher smoothing lag require
//     nSmoothingSymbolsVectors = 6;

	SNRs.push_back(3);
// 	SNRs.push_back(6);
//     SNRs.push_back(9);SNRs.push_back(12);SNRs.push_back(15);
// 	SNRs.push_back(18);SNRs.push_back(21);

    // BER and MSE computing
    symbolsDetectionWindowStart = trainSeqLength;
//     symbolsDetectionWindowStart = frameLength*3/10;    
    MSEwindowStart = 0;
    MSEwindowStart = frameLength*9/10;

    // results file name prefix
    sprintf(outputFileName,"res_");

    // alphabet is defined
    vector<vector<tBit> > alphabetBitSequences(2,vector<tBit>(1));
    alphabetBitSequences[0][0] = 0; alphabetBitSequences[1][0] = 1;
    vector<tSymbol> alphabetSymbols(2);
    alphabetSymbols[0] = -1; alphabetSymbols[1] = 1;
    alphabet = new Alphabet(alphabetBitSequences,alphabetSymbols);

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
    preamble = MatrixXd::Zero(1,1);
    preamble.resize(N,preambleLength);
    if(preamble.size()>0)
        preamble.setConstant(-1.0);
    
    // the frame length in bits is
    nBitsGenerated = (frameLength+nSmoothingSymbolsVectors)*alphabet->nBitsPerSymbol();
    
    // which symbols are to be taken into account when detecting
    isSymbolAccountedForDetection = vector<vector<bool> >(N,vector<bool>(frameLength));

    // the preamble symbols before symbolsDetectionWindowStart are ignored for detection
    for(int iTime=0;iTime<symbolsDetectionWindowStart;iTime++)
        for(int iInput=0;iInput<N;iInput++)
            isSymbolAccountedForDetection[iInput][iTime] = false;        
  
    for(int iTime=symbolsDetectionWindowStart;iTime<frameLength;iTime++)
        for(int iInput=0;iInput<N;iInput++)
            isSymbolAccountedForDetection[iInput][iTime] = true;   
    
#ifdef PRINT_SYMBOLS_ACCOUNTED_FOR_DETECTION
    cout << "isSymbolAccountedForDetection" << endl;
    Util::print(isSymbolAccountedForDetection);
#endif    
    
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
//     beforeRunStatUtilSeeds.reserve(nFrames);

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    presentFrameMSEtimeEvolution.resize(SNRs.size());
    MSEtimeEvolution.reserve(nFrames);
#endif

#ifdef KEEP_ALL_CHANNEL_MATRICES
	channelMatrices.reserve(nFrames);
#endif

#ifndef RANDOM_SEED
        // we don't want the same bits to be generated over and over
        randomGenerator.setSeed(0);
#endif

    channel = NULL;
    powerProfile = NULL;
    
#define XML_STRING_ATTRIBUTE(NAME,VALUE) "NAME=\"" << VALUE << "\"";    
    
    // data saving
    xmlFile.open("data.xml",ofstream::trunc);
    xmlFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    xmlFile << "<com>" << endl;
    xmlFile << "  <parameters>" << endl;
    xmlFile << "    <nInputs type=\"scalar\">" << N << "</nInputs>" << endl;
    xmlFile << "    <nOutputs type=\"scalar\">" << L << "</nOutputs>" << endl;
    xmlFile << "    <frameLength type=\"scalar\">" << frameLength << "</frameLength>" << endl;
    xmlFile << "    <channelOrder type=\"scalar\">" << m << "</channelOrder>" << endl;
    xmlFile << "    <smoothingLag type=\"scalar\">" << d << "</smoothingLag>" << endl;
    xmlFile << "    <trainingSeqLength type=\"scalar\">" << trainSeqLength << "</trainingSeqLength>" << endl;
    xmlFile << "    <preambleLength type=\"scalar\">" << preambleLength << "</preambleLength>" << endl;
    xmlFile << "  </parameters>" << endl;
}

BaseSystem::~BaseSystem()
{
    delete alphabet;
}

void BaseSystem::Simulate()
{

//     // for repeating simulations
//     randomGenerator.setSeed();
//     StatUtil::getRandomGenerator().setSeed();

    int iFrame = 0;
    while((iFrame<nFrames) && (!__done))
    {

        // the seeds are kept for saving later
        mainSeeds.push_back(randomGenerator.getSeed());
        statUtilSeeds.push_back(StatUtil::getRandomGenerator().getSeed());

        // bits are generated ...
        Bits bits(N,nBitsGenerated,randomGenerator);        

        // ... and then modulated by means of the alphabet
        MatrixXd symbolsWithoutPreamble = Modulator::modulate_eigen(bits,*alphabet);

        // the preamble is set before the symbols to be transmitted
        if(preamble.size()>0)
        {
            symbols.resize(preamble.rows(),preamble.cols()+symbolsWithoutPreamble.cols());
            symbols << preamble,symbolsWithoutPreamble;
        }else
            symbols = symbolsWithoutPreamble;

        // all the above symbols must be processed except those generated due to the smoothing
        iLastSymbolVectorToBeDetected = symbols.cols() - nSmoothingSymbolsVectors;

        BuildChannel();

#ifdef PRINT_PARAMETERS
        std::cout << "iLastSymbolVectorToBeDetected = " << iLastSymbolVectorToBeDetected << endl;
#endif

        // noise is generated according to the channel
//         noise = new NullNoise(L,channel->length());
//         noise = new ChannelDependentNoise(channel);
//         noise = new PowerProfileDependentNoise(L,channel->length(),*powerProfile);
		noise = new SingleUserPowerProfileDependentNoise(L,channel->length(),*powerProfile);

#ifdef EXPORT_REAL_DATA
            realSymbols = &symbols;
            realChannel = channel;
            realNoise = noise;
#endif

        for(iSNR=0;iSNR<SNRs.size();iSNR++)
        {
            cout << "SNR = " << SNRs[iSNR] << " [Frame " << iFrame << "]..." << endl;

            // noise SNR is set
            noise->setSNR(SNRs[iSNR],alphabet->variance());

#ifdef DEBUG
			cout << "noise is" << endl;
			noise->print();
			cout << endl;
#endif

            // transmission
            observations = channel->transmit(symbols,*noise);

             AddAlgorithms();

            // here the number of algoriths is known. So, the first iteration:
            if(iFrame==0 && iSNR==0)
                OnlyOnce();

            // algorithms are executed
            for(uint iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
            {
//                 // in order to repeat a concrete simulation...
//                 StatUtil::getRandomGenerator().setSeed();

                // the seed kept by the class StatUtil is saved
//                 presentFrameStatUtilSeeds(iSNR,iAlgorithm) = StatUtil::getRandomGenerator().getSeed();

                // if there is training sequence
                if(trainSeqLength!=0)
                    algorithms[iAlgorithm]->run(observations,noise->variances(),symbols.block(0,preambleLength,N,trainSeqLength));
                // if there is NOT training sequence
                else
                    algorithms[iAlgorithm]->run(observations,noise->variances());

                detectedSymbols = algorithms[iAlgorithm]->getDetectedSymbolVectors();
                
                pe = computeSER(symbols.block(0,preambleLength,N,frameLength),detectedSymbols,isSymbolAccountedForDetection);
				
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
		channel = NULL;
		
        delete noise;
		noise = NULL;
    } // while((iFrame<nFrames) && (!done))

    overallPeMatrix *= 1.0/iFrame;
    overallMseMatrix *= 1.0/iFrame;

    cout << "Overall SER:" << endl;
    Util::print(overallPeMatrix);

    cout << "Overall MSE:" << endl;
    Util::print(overallMseMatrix);
    
    // data saving
    xmlFile << "</com>" << endl;
    xmlFile.close();
}

void BaseSystem::OnlyOnce()
{
    overallPeMatrix = MatrixXd::Zero(SNRs.size(),algorithms.size());
    presentFramePe = MatrixXd::Zero(SNRs.size(),algorithms.size());

    overallMseMatrix = MatrixXd::Zero(SNRs.size(),algorithms.size());
    presentFrameMSE = MatrixXd::Zero(SNRs.size(),algorithms.size());

    // Pe evolution
    for(uint i=0;i<SNRs.size();i++)
    {
        overallPeTimeEvolution[i] = MatrixXd(algorithms.size(),frameLength);
        overallErrorsNumberTimeEvolution[i] = MatrixXi::Zero(algorithms.size(),frameLength);
    }

    // we fill the vector with the names of the algorithms
    for(uint iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
        algorithmsNames.push_back(algorithms[iAlgorithm]->getName());

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    for(uint i=0;i<SNRs.size();i++)
        presentFrameMSEtimeEvolution[i] = MatrixXd(algorithms.size(),frameLength);
#endif

//     presentFrameStatUtilSeeds = LaGenMatLongInt(SNRs.size(),algorithms.size());
}

void BaseSystem::BeforeEndingAlgorithm(int iAlgorithm)
{
    mse = algorithms[iAlgorithm]->MSE(channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1));

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    VectorXd mseAlongTime = TransmissionUtil::MSEalongTime(algorithms[iAlgorithm]->getEstimatedChannelMatrices(),0,frameLength-1,channel->range_eigen(preambleLength,preambleLength+frameLength-1),0,frameLength-1);
    for(int ik=0;ik<frameLength;ik++)
        presentFrameMSEtimeEvolution[iSNR](iAlgorithm,ik) = mseAlongTime(ik);
#endif

    cout << algorithms[iAlgorithm]->getName() << ": Pe = " << pe << " , MSE = " << mse << endl;

    // the error probability is accumulated
    overallPeMatrix(iSNR,iAlgorithm) += pe;
    presentFramePe(iSNR,iAlgorithm) = pe;

    // and the MSE
    overallMseMatrix(iSNR,iAlgorithm) += mse;
    presentFrameMSE(iSNR,iAlgorithm) = mse;

    // Pe evolution
    MatrixXd transmittedSymbols = symbols.block(0,preambleLength,N,frameLength);

    if(detectedSymbols.rows()!=0)
    {
        for(int k=0;k<frameLength;k++)
            for(int iUser=0;iUser<N;iUser++)
                if(detectedSymbols(iUser,k)!=transmittedSymbols(iUser,k))
                    overallErrorsNumberTimeEvolution[iSNR](iAlgorithm,k)++;
    }
}

void BaseSystem::BeforeEndingFrame(int iFrame)
{
    // pe
    peMatrices.push_back(presentFramePe);
    Util::matricesVectorToOctaveFileStream(peMatrices,"pe",f);

    // MSE
    MSEMatrices.push_back(presentFrameMSE);
    Util::matricesVectorToOctaveFileStream(MSEMatrices,"mse",f);

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    MSEtimeEvolution.push_back(presentFrameMSEtimeEvolution);
    Util::matricesVectorsVectorToOctaveFileStream(MSEtimeEvolution,"MSEtimeEvolution",f);
#endif

#ifdef KEEP_ALL_CHANNEL_MATRICES
	channelMatrices.push_back(channel->range(preambleLength,iLastSymbolVectorToBeDetected));
#endif

    // seeds just before the run of the algorithms
//     beforeRunStatUtilSeeds.push_back(presentFrameStatUtilSeeds);
//     Util::matricesVectorToOctaveFileStream(beforeRunStatUtilSeeds,"beforeRunStatUtilSeeds",f);

//     for(uint iSNR=0;iSNR<SNRs.size();iSNR++)
//         for(uint i=0;i<algorithmsNames.size();i++)
//             for(int j=0;j<frameLength;j++)
//                 overallPeTimeEvolution[iSNR](i,j) = (double) overallErrorsNumberTimeEvolution[iSNR](i,j) / (double) (N*(iFrame+1));
//     Util::matricesVectorToOctaveFileStream(overallPeTimeEvolution,"peTimeEvolution",f);

    Util::scalarToOctaveFileStream(iFrame+1,"nFrames",f);

    Util::stringsVectorToOctaveFileStream(algorithmsNames,"algorithmsNames",f);
    Util::scalarToOctaveFileStream(L,"L",f);
    Util::scalarToOctaveFileStream(N,"N",f);
    Util::scalarToOctaveFileStream(m,"m",f);
    Util::scalarToOctaveFileStream(frameLength,"frameLength",f);
    Util::scalarToOctaveFileStream(trainSeqLength,"trainSeqLength",f);
    Util::scalarToOctaveFileStream(d,"d",f);
    Util::scalarToOctaveFileStream(symbolsDetectionWindowStart,"symbolsDetectionWindowStart",f);
    Util::scalarToOctaveFileStream(MSEwindowStart,"MSEwindowStart",f);
    Util::scalarsVectorToOctaveFileStream(SNRs,"SNRs",f);
    Util::matrixToOctaveFileStream(preamble,"preamble",f);
    Util::scalarToOctaveFileStream(nSmoothingSymbolsVectors,"nSmoothingSymbolsVectors",f);    
    Util::scalarToOctaveFileStream(preambleLength,"preambleLength",f);
    Util::scalarsVectorToOctaveFileStream(mainSeeds,"mainSeeds",f);
    Util::scalarsVectorToOctaveFileStream(statUtilSeeds,"statUtilSeeds",f);
    Util::matricesVectorToOctaveFileStream(channel->range(preambleLength,iLastSymbolVectorToBeDetected),"channel",f);
	
#ifdef KEEP_ALL_CHANNEL_MATRICES
	Util::matricesVectorsVectorToOctaveFileStream(channelMatrices,"channels",f);
#endif

    Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*channel).name())),"channelClass",f);
    Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*noise).name())),"noiseClass",f);
    Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*this).name())),"systemClass",f);

    if(powerProfile!=NULL)
    {
        Util::scalarsVectorToOctaveFileStream(powerProfile->tapsPowers(),"powerProfileVariances",f);
        Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*powerProfile).name())),"powerProfileClass",f);
    }
    
}

double BaseSystem::computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask)
{
    if(detectedSymbols.rows() == 0)
        return 0.0;

    if(sourceSymbols.rows()!= detectedSymbols.rows() || detectedSymbols.rows()!= mask.size())
    {
        cout << "sourceSymbols.rows() = " << sourceSymbols.rows() << " detectedSymbols.rows() = " << detectedSymbols.rows() << " mask.size() = " << mask.size() << endl;
        throw RuntimeException("TransmissionUtil::computeSER: matrix row numbers differ.");
    }

    if(sourceSymbols.cols()!= detectedSymbols.cols() || detectedSymbols.cols()!= mask[0].size())
    {
        cout << "sourceSymbols.cols() = " << sourceSymbols.cols() << " detectedSymbols.cols() = " << detectedSymbols.cols() << " mask.size() = " << mask.size() << endl;    
      throw RuntimeException("TransmissionUtil::computeSER: matrix column numbers differ.");
    }
        
#ifdef PRINT_INFO
    cout << "source symbols" << endl << sourceSymbols << "detected symbols" << endl << detectedSymbols << "mask" << endl;
    Util::print(mask);
#endif

    uint iBestPermutation = 0;
    vector<int> bestPermutationSigns(sourceSymbols.rows(),1);

    // max number of errors
    int minErrors = sourceSymbols.rows()*sourceSymbols.cols()*alphabet->nBitsPerSymbol();
    
    uint nAccountedSymbols = 0;
    uint iInput;

    for(uint iPermut=0;iPermut<permutations.size();iPermut++)
    {
        int permutationErrors = 0;
        
        for(uint iStream=0;iStream<permutations[iPermut].size();iStream++)
        {
            iInput = permutations[iPermut][iStream];
          
            int errorsInverting=0,errorsWithoutInverting=0;
            
            for(uint iTime=0;iTime<static_cast<uint> (sourceSymbols.cols());iTime++)
            {
                // if this symbol is not accounted for
                if(!mask[iStream][iTime])
                    continue;

                // if the symbols differ, an error happened...
                errorsWithoutInverting += sourceSymbols(iStream,iTime) != detectedSymbols(iInput,iTime);
                
                // ...unless there the symbol sign needs to be switched because of the ambiguity
                errorsInverting += sourceSymbols(iStream,iTime) != alphabet->opposite(detectedSymbols(iInput,iTime));
                
                nAccountedSymbols++;
            }              
            
            if(errorsWithoutInverting<errorsInverting)
            {
                permutationErrors += errorsWithoutInverting;
                bestPermutationSigns[iStream] = 1;
            }
            else
            {
                permutationErrors += errorsInverting;
                bestPermutationSigns[iStream] = -1;
            }
        } // for(uint iStream=0;iStream<permutations[iPermut].size();iStream++)
        
        if(permutationErrors<minErrors)
        {
            minErrors = permutationErrors;
            iBestPermutation = iPermut;
        }
    }
    
    nAccountedSymbols /= permutations.size();
    
    return (double)minErrors/(double)(nAccountedSymbols);
}