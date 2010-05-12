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
#include <bashcolors.h>
#include <typeinfo>
#include <string.h>
#include <assert.h>

#include <SingleUserPowerProfileDependentNoise.h>

#define DATE_LENGTH 100

#define EXPORT_REAL_DATA

// #define PRINT_PARAMETERS
// #define PRINT_SYMBOLS_ACCOUNTED_FOR_DETECTION
// #define PRINT_SYMBOLS_ACCOUNTED_FOR_DETECTION_PER_FRAME

// #define PRINT_COMPUTE_SER_INFO
// #define PRINT_BEST_PERMUATION_WHEN_COMPUTING_SER
// #define PRINT_BEST_PERMUATION_ERRORS

// #define STOP_AFTER_EACH_FRAME
// #define STOP_AFTER_EACH_SNR

#define SAVE_SEEDS
#define LOAD_SEEDS

// #define DEBUG
// #define DEBUG2
// #define DEBUG_SER_WITHOUT_SOLVING_AMBIGUITY

using namespace std;

#ifdef EXPORT_REAL_DATA
    MIMOChannel *realChannel;
    MatrixXd *realSymbols;
    Noise *realNoise;
#endif

BaseSystem::BaseSystem()
{
    // GLOBAL PARAMETERS

// ------------------------ iswcs 2010 ----------------------
//     nFrames = 100;
//     L=3,N=3,frameLength=300;
//     m = 3;
//     d = m - 1;
//     trainSeqLength = 10;
//     preambleLength = 10;
//   
//     // the algorithms with the higher smoothing lag require
//     nSmoothingSymbolsVectors = 10;

// ---------------------------- tesis ------------------------

//     nFrames = 10;
//     L=3,N=2,frameLength=300;
//     m = 3;
//     d = m - 1;
//     trainSeqLength = 10;
//     preambleLength = 10;
//   
//     // the algorithms with the higher smoothing lag require
//     nSmoothingSymbolsVectors = 10;

// --------------------------- CDMA -------------------------

	nFrames = 10000;
// 	nFrames = 1;
// 	nFrames = 200;
//     L=3,N=2,frameLength=300;
//     L=8,N=3,frameLength=300;
    L=8,N=3,frameLength=1000;
    m = 1;
    d = m - 1;
    trainSeqLength = 0;
    preambleLength = 0;
    
    // the algorithms with the higher smoothing lag require
    nSmoothingSymbolsVectors = 6;

//   SNRs.push_back(0);
  SNRs.push_back(3);
  SNRs.push_back(6);
  SNRs.push_back(9);SNRs.push_back(12);SNRs.push_back(15);
//   SNRs.push_back(18);SNRs.push_back(21);

    // BER and MSE computing
    symbolsDetectionWindowStart = trainSeqLength;
//     symbolsDetectionWindowStart = frameLength*3/10; 

    MSEwindowStart = frameLength*9/10;
//     MSEwindowStart = 0;

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
    permutations = Util::permutations(firstPermutation,N);
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

#ifdef KEEP_ALL_CHANNEL_ESTIMATIONS
	channelEstimations.reserve(nFrames);
#endif

#ifndef RANDOM_SEED
        // we don't want the same bits to be generated over and over
        randomGenerator.setSeed(0);
#endif

    channel = NULL;
    powerProfile = NULL;
    
// #define XML_STRING_ATTRIBUTE(NAME,VALUE) "NAME=\"" << VALUE << "\"";    
//     
//     // data parameters saving
//     xmlFile.open("data.xml",ofstream::trunc);
//     xmlFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
//     xmlFile << "<com>" << endl;
//     xmlFile << "  <parameters>" << endl;
//     xmlFile << "    <nInputs type=\"scalar\">" << N << "</nInputs>" << endl;
//     xmlFile << "    <nOutputs type=\"scalar\">" << L << "</nOutputs>" << endl;
//     xmlFile << "    <frameLength type=\"scalar\">" << frameLength << "</frameLength>" << endl;
//     xmlFile << "    <channelOrder type=\"scalar\">" << m << "</channelOrder>" << endl;
//     xmlFile << "    <smoothingLag type=\"scalar\">" << d << "</smoothingLag>" << endl;
//     xmlFile << "    <trainingSeqLength type=\"scalar\">" << trainSeqLength << "</trainingSeqLength>" << endl;
//     xmlFile << "    <preambleLength type=\"scalar\">" << preambleLength << "</preambleLength>" << endl;
//     xmlFile << "  </parameters>" << endl;
}

BaseSystem::~BaseSystem()
{
    delete alphabet;
}

void BaseSystem::Simulate()
{

#ifdef LOAD_SEEDS
    // for repeating simulations
    randomGenerator.setSeed(1251453452);
    StatUtil::getRandomGenerator().setSeed(4117327754);

//     randomGenerator.setSeed(3013783408);
//     StatUtil::getRandomGenerator().setSeed(1878942058);

    cout << COLOR_LIGHT_BLUE << "seeds are being loaded..." << COLOR_NORMAL << endl;
#endif

    iFrame = 0;
    while((iFrame<nFrames) && (!__done))
    {

#ifdef SAVE_SEEDS
        // the seeds are kept for saving later
        mainSeeds.push_back(randomGenerator.getSeed());
        statUtilSeeds.push_back(StatUtil::getRandomGenerator().getSeed());
#endif

        // bits are generated ...
        Bits generatedBits(N,nBitsGenerated,randomGenerator);

		// differential modulation
// 		bits = bits.differentialEncoding();

        // ... and then modulated by means of the alphabet
        MatrixXd symbolsWithoutPreamble = Modulator::modulate(generatedBits,*alphabet);

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

#ifdef PRINT_SYMBOLS_ACCOUNTED_FOR_DETECTION_PER_FRAME
    cout << "isSymbolAccountedForDetection" << endl;
    Util::print(isSymbolAccountedForDetection);
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
            cout << COLOR_FRAME_NUMBER_SNR << "SNR = " << COLOR_NORMAL << SNRs[iSNR] << COLOR_FRAME_NUMBER_SNR << " [Frame " << COLOR_NORMAL << iFrame << COLOR_FRAME_NUMBER_SNR << "]..." << COLOR_NORMAL << endl;

            // noise SNR is set
            noise->setSNR(SNRs[iSNR],alphabet->variance());

#ifdef PRINT_NOISE
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
            for(iAlgorithm=0;iAlgorithm<algorithms.size();iAlgorithm++)
            {
                // if there is training sequence
                if(trainSeqLength!=0)
                    algorithms[iAlgorithm]->run(observations,noise->variances(),symbols.block(0,preambleLength,N,trainSeqLength));
                // if there is NOT training sequence
                else
                    algorithms[iAlgorithm]->run(observations,noise->variances());

                detectedSymbols = algorithms[iAlgorithm]->getDetectedSymbolVectors();

                pe = computeSER(symbols.block(0,preambleLength,N,frameLength),detectedSymbols,isSymbolAccountedForDetection,_iBestPermutation,_bestPermutationSigns);

                BeforeEndingAlgorithm();

                delete algorithms[iAlgorithm];
            }

            algorithms.clear();
			
#ifdef STOP_AFTER_EACH_SNR
		getchar();
#endif
        } // for(int iSNR=0;iSNR<SNRs.size();iSNR++)

        f.open(outputFileName,ofstream::trunc);

        BeforeEndingFrame();
        
		f.close();

        // ---------------------------------------------------------

        iFrame++;

        delete channel;
		channel = NULL;
		
        delete noise;
		noise = NULL;
		
#ifdef STOP_AFTER_EACH_FRAME
		getchar();
#endif
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

#ifdef KEEP_ALL_CHANNEL_ESTIMATIONS
  // channel estimations
  presentFrameChannelMatrixEstimations = std::vector<std::vector<std::vector<MatrixXd> > >(SNRs.size(),std::vector<std::vector<MatrixXd> >(algorithms.size()));
#endif
}

void BaseSystem::BeforeEndingAlgorithm()
{
//     mse = algorithms[iAlgorithm]->MSE(channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1),permutations[_iBestPermutation],_bestPermutationSigns);
//     mse = algorithms[iAlgorithm]->MSE(channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1));
	
	// the channel matrices estimated by the algorithm are stored
	std::vector<MatrixXd> estimatedChannelMatrices = algorithms[iAlgorithm]->getEstimatedChannelMatrices();
	
	if(estimatedChannelMatrices.size()!=0)
	{
	  std::vector<MatrixXd> toCheckEstimatedChannelMatrices(estimatedChannelMatrices.begin()+MSEwindowStart,estimatedChannelMatrices.end());
// 	  cout << "estimatedChannelMatrices.size() = " << estimatedChannelMatrices.size() << " toCheckEstimatedChannelMatrices.size() = " << toCheckEstimatedChannelMatrices.size() << endl;
// 	  cout << "channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1).size() = " << channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1).size() << endl;

	  mse = computeMSE(channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1),toCheckEstimatedChannelMatrices);
	}else
	{
	  // it should be zero
	  mse = computeMSE(channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1),estimatedChannelMatrices);
	}
	
#ifdef MSE_TIME_EVOLUTION_COMPUTING
    VectorXd mseAlongTime = TransmissionUtil::MSEalongTime(algorithms[iAlgorithm]->getEstimatedChannelMatrices(),0,frameLength-1,channel->range(preambleLength,preambleLength+frameLength-1),0,frameLength-1);
    for(int ik=0;ik<frameLength;ik++)
        presentFrameMSEtimeEvolution[iSNR](iAlgorithm,ik) = mseAlongTime(ik);
#endif

#ifdef KEEP_ALL_CHANNEL_ESTIMATIONS
  // we get the channel matrices estimated by this algorithm
  vector<MatrixXd> thisAlgorithmEstimatedChannelMatrices = algorithms[iAlgorithm]->getEstimatedChannelMatrices();

  // if none, that meaning the algorithm does not performa channel matrix estimation,...
  if(thisAlgorithmEstimatedChannelMatrices.size()==0)
	// we generate a sequence of matrices
	thisAlgorithmEstimatedChannelMatrices = vector<MatrixXd>(iLastSymbolVectorToBeDetected-preambleLength,MatrixXd::Zero(channel->channelCoefficientsMatrixRows(),channel->channelCoefficientsMatrixCols()));
  else
	if(thisAlgorithmEstimatedChannelMatrices.size()!=iLastSymbolVectorToBeDetected-preambleLength)
	  throw RuntimeException("BaseSystem::BeforeEndingAlgorithm: the number of channel matrices estimated by the algorithm is not the expected.");

  presentFrameChannelMatrixEstimations[iSNR][iAlgorithm] = thisAlgorithmEstimatedChannelMatrices;
#endif

    cout << COLOR_GREEN << algorithms[iAlgorithm]->getName() << COLOR_NORMAL << ": Pe = " << pe << " , MSE = " << mse << endl;

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

void BaseSystem::BeforeEndingFrame()
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
	channelMatrices.push_back(channel->range(preambleLength,iLastSymbolVectorToBeDetected-1));
#endif

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
	
#ifdef KEEP_ALL_CHANNEL_MATRICES
	Util::matricesVectorsVectorToOctaveFileStream(channelMatrices,"channels",f);
#else
    Util::matricesVectorToOctaveFileStream(channel->range(preambleLength,iLastSymbolVectorToBeDetected),"channel",f);
#endif

    Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*channel).name())),"channelClass",f);
    Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*noise).name())),"noiseClass",f);
    Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*this).name())),"systemClass",f);

    if(powerProfile!=NULL)
    {
        Util::scalarsVectorToOctaveFileStream(powerProfile->tapsPowers(),"powerProfileVariances",f);
        Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*powerProfile).name())),"powerProfileClass",f);
    }
    
#ifdef KEEP_ALL_CHANNEL_ESTIMATIONS
	channelEstimations.push_back(presentFrameChannelMatrixEstimations);
	Util::matricesVectorsVectorsVectoresVectorToOctaveFileStream(channelEstimations,"channelEstimations",f);
#endif
}

double BaseSystem::computeSER(const MatrixXd& sourceSymbols, const MatrixXd& detectedSymbols, const std::vector< std::vector< bool > >& mask, uint &iBestPermutation, std::vector< int > &bestPermutationSigns)
{
	// best permutation and its signs are initialized...
    iBestPermutation = 0;
	bestPermutationSigns = vector<int>(sourceSymbols.rows(),1);
	
	//...even if there are no symbols detected
    if(detectedSymbols.rows() == 0)
        return 0.0;

    if(sourceSymbols.rows()!= detectedSymbols.rows() || detectedSymbols.rows()!= mask.size())
    {
        cout << "sourceSymbols.rows() = " << sourceSymbols.rows() << " detectedSymbols.rows() = " << detectedSymbols.rows() << " mask.size() = " << mask.size() << endl;
        throw RuntimeException("BaseSystem::computeSER: matrix row numbers differ.");
    }

    if(sourceSymbols.cols()!= detectedSymbols.cols() || detectedSymbols.cols()!= mask[0].size())
    {
        cout << "sourceSymbols.cols() = " << sourceSymbols.cols() << " detectedSymbols.cols() = " << detectedSymbols.cols() << " mask.size() = " << mask.size() << endl;    
      throw RuntimeException("BaseSystem::computeSER: matrix column numbers differ.");
    }
        
#ifdef PRINT_COMPUTE_SER_INFO
    cout << "source symbols" << endl << sourceSymbols << endl << endl << "detected symbols" << endl << detectedSymbols << endl << endl << "mask" << endl;
    Util::print(mask);
#endif

//     iBestPermutation = 0;
    vector<int> thisPermutationSigns(sourceSymbols.rows());
// 	bestPermutationSigns = vector<int>(sourceSymbols.rows());

    // max number of errors
    int minErrors = sourceSymbols.rows()*sourceSymbols.cols();
    
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
                errorsWithoutInverting += (sourceSymbols(iStream,iTime) != detectedSymbols(iInput,iTime));
                
                // ...unless the symbol sign needs to be switched because of the ambiguity
                errorsInverting += (sourceSymbols(iStream,iTime) != alphabet->opposite(detectedSymbols(iInput,iTime)));
                
                nAccountedSymbols++;
            }              

#ifdef PRINT_BEST_PERMUATION_ERRORS
			if(iPermut==0)
			{
			cout << "====================" << endl;
			cout << "iStream = " << iStream << endl;
			cout << "errorsWithoutInverting = " << errorsWithoutInverting << endl;
			cout << "errorsInverting = " << errorsInverting << endl;
			cout << "====================" << endl;
			}
#endif

            if(errorsWithoutInverting<errorsInverting)
            {
                permutationErrors += errorsWithoutInverting;
                thisPermutationSigns[iStream] = 1;
            }
            else
            {
                permutationErrors += errorsInverting;
                thisPermutationSigns[iStream] = -1;
            }
        } // for(uint iStream=0;iStream<permutations[iPermut].size();iStream++)
        
        if(permutationErrors<minErrors)
        {
            minErrors = permutationErrors;
            iBestPermutation = iPermut;
			bestPermutationSigns = thisPermutationSigns;
        }
    }
    
#ifdef PRINT_BEST_PERMUATION_WHEN_COMPUTING_SER
	cout << "best permutation is " << iBestPermutation << endl;
	cout << "its signs" << endl;
	Util::print(bestPermutationSigns);
	cout << endl;
#endif
    
	assert(nAccountedSymbols % permutations.size() == 0);
    
    nAccountedSymbols /= permutations.size();
    
    // if all the symbols were masked
    if(nAccountedSymbols==0)
      return 0.0;
    else
      return (double)minErrors/(double)(nAccountedSymbols);
}

double BaseSystem::computeMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices) const
{
    int nRealChannelMatrices = realChannelMatrices.size();
    int nEstimatedChannelMatrices = estimatedChannelMatrices.size();

    // if the algorithm didn't perform channel estimation
    if(nEstimatedChannelMatrices==0)
        return -1.0;

    if(nRealChannelMatrices!=nEstimatedChannelMatrices)
        throw RuntimeException("BaseSystem::computeMSE: number of real channel matrices doesn't match that of the estimated.");

    double mse = 0;

	
// 	cout << "realChannelMatrices.at(0) = " << endl << realChannelMatrices.at(0) << endl;
// 	cout << "estimatedChannelMatrices.at(0) = " << endl << estimatedChannelMatrices.at(0) << endl;
	
	for(int i=0;i<nRealChannelMatrices;i++)
	{
		// the square error committed by the estimated matrix is normalized by the squared Frobenius norm
		// (i.e. the sum of all the elements squared) of the real channel matrix
		// also notice that if the channel is Sparkling memory, the channel matrices of the real channel may have different sizes
		mse += Util::squareErrorPaddingWithZeros(realChannelMatrices.at(i),estimatedChannelMatrices.at(i))/pow(realChannelMatrices.at(i).norm(),2.0);
	}

    return mse/(double)nRealChannelMatrices;
}

double BaseSystem::computeMSE(const vector<MatrixXd> &realchannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices,const vector<uint> &bestPermutation,const vector<int> &bestPermutationSigns) const
{
  vector<uint> realChannelMatricesPermutation = Util::computeInversePermutation(bestPermutation);

  // the signs permutation is given permuted WITH RESPECT TO THE ESTIMATED CHANNEL MATRICES
  // we have to permute them according to the best permutation with respect to the real channel matrices
  vector<int> realChannelMatricesSignsPermutation = Util::applyPermutation(Util::applyPermutation(bestPermutationSigns,realChannelMatricesPermutation),realChannelMatricesPermutation);

  vector<MatrixXd> permutedRealChannelMatrices(realchannelMatrices.size());
  
  for(uint i=0;i<realchannelMatrices.size();i++)
  {
	permutedRealChannelMatrices[i] = Util::applyPermutationOnColumns(realchannelMatrices[i],realChannelMatricesPermutation,realChannelMatricesSignsPermutation);
// 	cout << "realchannelMatrices[i] = " << endl << realchannelMatrices[i] << endl;
// 	cout << "permutedRealChannelMatrices[i] = " << endl << permutedRealChannelMatrices[i] << endl;
  }

  return BaseSystem::computeMSE(permutedRealChannelMatrices,estimatedChannelMatrices);
}

double BaseSystem::computeSERwithoutSolvingAmbiguity(const MatrixXd& sourceSymbols, const MatrixXd& detectedSymbols, const std::vector< std::vector< bool > >& mask) const
{
#ifdef DEBUG
  cout << "sourceSymbols.cols() = " << sourceSymbols.cols() << " detectedSymbols.cols() = " << detectedSymbols.cols() << " mask[0].size() = " << mask[0].size() << endl;
#endif
  
  if(detectedSymbols.rows() == 0)
	  return -1.0;

  if(sourceSymbols.rows()!= detectedSymbols.rows() || detectedSymbols.rows()!= mask.size())
  {
	  cout << "sourceSymbols.rows() = " << sourceSymbols.rows() << " detectedSymbols.rows() = " << detectedSymbols.rows() << " mask.size() = " << mask.size() << endl;
	  throw RuntimeException("BaseSystem::computeSERwithoutSolvingAmbiguity: matrix row numbers differ.");
  }

  if(sourceSymbols.cols()!= detectedSymbols.cols() || detectedSymbols.cols()!= mask[0].size())
  {
	  cout << "sourceSymbols.cols() = " << sourceSymbols.cols() << " detectedSymbols.cols() = " << detectedSymbols.cols() << " mask.size() = " << mask.size() << endl; 
	throw RuntimeException("BaseSystem::computeSERwithoutSolvingAmbiguity: matrix column numbers differ.");
  }

#ifdef DEBUG_SER_WITHOUT_SOLVING_AMBIGUITY
  if(iAlgorithm==3 && iSNR==4)
  {
	cout << "sourceSymbols" << endl << sourceSymbols << endl;
	cout << "detectedSymbols" << endl << detectedSymbols << endl;
	Util::print(mask);
	cout << endl;
  }
#endif

  // max number of errors
  uint errors = 0;

  uint nAccountedSymbols = 0;

  for(uint iStream=0;iStream<N;iStream++)
  {
	  for(uint iTime=0;iTime<static_cast<uint> (sourceSymbols.cols());iTime++)
	  {
		  // if this symbol is not accounted for
		  if(!mask[iStream][iTime])
			  continue;

		  // if the symbols differ, an error happened...
		  errors += (sourceSymbols(iStream,iTime) != detectedSymbols(iStream,iTime));
		  
		  nAccountedSymbols++;
	  }              
  } // for(uint iStream=0;iStream<permutations[iPermut].size();iStream++)

  // if all the symbols were masked
  if(nAccountedSymbols==0)
	return 0.0;
  else
	return (double)errors/(double)(nAccountedSymbols);
}