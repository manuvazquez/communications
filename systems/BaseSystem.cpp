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
#include <stdlib.h>

#include <SingleUserPowerProfileDependentNoise.h>

extern uint32_t __mainSeedPassed;
extern uint32_t __statUtilSeedPassed;
extern uint __nFramesPassed;
extern bool __randomSeedHasBeenPassed;
extern bool __nFramesHasBeenPassed;

#define DATE_LENGTH 100

// #define EXPORT_REAL_DATA
#define EXPORT_REAL_CHANNEL_ORDER

// #define PRINT_PARAMETERS
// #define PRINT_SYMBOLS_ACCOUNTED_FOR_DETECTION
// #define PRINT_SYMBOLS_ACCOUNTED_FOR_DETECTION_PER_FRAME

// #define PRINT_COMPUTE_SER_INFO
// #define PRINT_BEST_PERMUATION_WHEN_COMPUTING_SER
// #define PRINT_BEST_PERMUATION_ERRORS

// #define STOP_AFTER_EACH_FRAME
// #define STOP_AFTER_EACH_SNR

#define SAVE_SEEDS
// #define LOAD_SEEDS

// #define DEBUG
// #define DEBUG2
// #define DEBUG_SER_WITHOUT_SOLVING_AMBIGUITY

#ifdef EXPORT_REAL_DATA
    MIMOChannel *realChannel;
    MatrixXd *realSymbols;
    Noise *realNoise;
#endif

#ifdef EXPORT_REAL_CHANNEL_ORDER
	int realChannelOrder;
#endif

// #define TYPEID(x) strpbrk(typeid(x).name(),"_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")

BaseSystem::BaseSystem()
{
    // GLOBAL PARAMETERS

// ------------------------ iswcs 2010 ----------------------
	_nFrames = 500;
	_L=3,_N=3,_frameLength=300;
	_m = 4;
	_d = _m - 1;

	_trainSeqLength = 15;

	_preambleLength = 10;

	// the algorithms with the higher smoothing lag require
	_nSmoothingSymbolsVectors = 10;

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

// // 	_nFrames = 10000;
// 	_nFrames = 500;
// 	_L=8,_N=3,_frameLength=1000;
// 	_m = 1;
// 	_d = _m - 1;
// 	_trainSeqLength = 0;
// 	_preambleLength = 0;
// 
// 	// the algorithms with the higher smoothing lag require
// 	_nSmoothingSymbolsVectors = 6;

	_SNRs.push_back(0);
	_SNRs.push_back(3);
	_SNRs.push_back(6);
	_SNRs.push_back(9);
	_SNRs.push_back(12);
	_SNRs.push_back(15);
// 	_SNRs.push_back(18);
// 	_SNRs.push_back(21);
// 	_SNRs.push_back(24);
// 	_SNRs.push_back(27);

    // BER and MSE computing
    _symbolsDetectionWindowStart = _trainSeqLength;
//     symbolsDetectionWindowStart = frameLength*3/10; 

    _MSEwindowStart = _frameLength*9/10;
//     MSEwindowStart = 0;

	// ---------------------------------------------------------------------------------------------

#ifdef EXPORT_REAL_CHANNEL_ORDER
	realChannelOrder = _m;
#endif
    // results file name prefix
    sprintf(_outputFileName,"res_");

    // alphabet is defined
    vector<vector<tBit> > alphabetBitSequences(2,vector<tBit>(1));
    alphabetBitSequences[0][0] = 0; alphabetBitSequences[1][0] = 1;
    vector<tSymbol> alphabetSymbols(2);
    alphabetSymbols[0] = -1; alphabetSymbols[1] = 1;
    _alphabet = new Alphabet(alphabetBitSequences,alphabetSymbols);

    // host name is concatenated into the file name
    char hostname[HOSTNAME_LENGTH];
    gethostname(hostname,HOSTNAME_LENGTH);
    strcat(_outputFileName,hostname);

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
    strcat(_outputFileName,"_");
    strcat(_outputFileName,presentTimeString);
	
	// the name for the temporal file is obtained from the final one
	 strcpy(_tempOutputFileName,_outputFileName);
	 strcat(_tempOutputFileName,"_tmp");

    // a specific preamble is generated...
    _preamble = MatrixXd::Zero(1,1);
    _preamble.resize(_N,_preambleLength);
    if(_preamble.size()>0)
        _preamble.setConstant(-1.0);
    
    // the frame length in bits is
    _nBitsGenerated = (_frameLength+_nSmoothingSymbolsVectors)*_alphabet->nBitsPerSymbol();
    
    // which symbols are to be taken into account when detecting
    _isSymbolAccountedForDetection = vector<vector<bool> >(_N,vector<bool>(_frameLength));

    // the preamble symbols before symbolsDetectionWindowStart are ignored for detection
    for(int iTime=0;iTime<_symbolsDetectionWindowStart;iTime++)
        for(uint iInput=0;iInput<_N;iInput++)
            _isSymbolAccountedForDetection[iInput][iTime] = false;        
  
    for(int iTime=_symbolsDetectionWindowStart;iTime<_frameLength;iTime++)
        for(uint iInput=0;iInput<_N;iInput++)
            _isSymbolAccountedForDetection[iInput][iTime] = true;   
    
#ifdef PRINT_SYMBOLS_ACCOUNTED_FOR_DETECTION
    cout << "isSymbolAccountedForDetection" << endl;
    Util::print(isSymbolAccountedForDetection);
#endif    
    
    // ambiguity resolution
    uint *firstPermutation = new uint[_N];
    for(uint i=0;i<_N;i++) firstPermutation[i] = i;
    _permutations = Util::permutations(firstPermutation,_N);
    delete[] firstPermutation;

    _peMatrices.reserve(_nFrames);
    _MSEMatrices.reserve(_nFrames);

    _overallPeTimeEvolution.resize(_SNRs.size());
    _overallErrorsNumberTimeEvolution.resize(_SNRs.size());

    _mainSeeds.reserve(_nFrames);
    _statUtilSeeds.reserve(_nFrames);
//     beforeRunStatUtilSeeds.reserve(nFrames);

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    presentFrameMSEtimeEvolution.resize(SNRs.size());
    MSEtimeEvolution.reserve(nFrames);
#endif

#ifdef KEEP_ALL_CHANNEL_MATRICES
	_channelMatrices.reserve(_nFrames);
#endif

#ifdef KEEP_ALL_CHANNEL_ESTIMATIONS
	channelEstimations.reserve(_nFrames);
#endif

#ifndef RANDOM_SEED
        // we don't want the same bits to be generated over and over
        _randomGenerator.setSeed(0);
#endif

    _channel = NULL;
    _powerProfile = NULL;
    
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
    delete _alphabet;
}

void BaseSystem::simulate()
{

// if a number of frames has been passed...
if(__nFramesHasBeenPassed)
{
	_nFrames = __nFramesPassed;
	cout << COLOR_LIGHT_BLUE << _nFrames << " frames are going to be simulated." << COLOR_NORMAL << endl;
}

#ifdef LOAD_SEEDS
    // for repeating simulations
	
	if(__randomSeedHasBeenPassed)
	{
		_randomGenerator.setSeed(__mainSeedPassed);
		StatUtil::getRandomGenerator().setSeed(__statUtilSeedPassed);
	}else
	{
		_randomGenerator.setSeed(186951016);
		StatUtil::getRandomGenerator().setSeed(3731316448);
	}

    cout << COLOR_LIGHT_BLUE << "seeds are being loaded..." << COLOR_NORMAL << endl;
	cout << COLOR_LIGHT_BLUE << "\t " << _randomGenerator.getSeed() << endl << "\t " << StatUtil::getRandomGenerator().getSeed() << COLOR_NORMAL << endl;
#endif

    _iFrame = 0;
    while((_iFrame<_nFrames) && (!__done))
    {

#ifdef SAVE_SEEDS
        // the seeds are kept for saving later
        _mainSeeds.push_back(_randomGenerator.getSeed());
        _statUtilSeeds.push_back(StatUtil::getRandomGenerator().getSeed());
#endif

        // bits are generated ...
        Bits generatedBits(_N,_nBitsGenerated,_randomGenerator);

		// differential modulation
// 		bits = bits.differentialEncoding();

        // ... and then modulated by means of the alphabet
        MatrixXd symbolsWithoutPreamble = Modulator::modulate(generatedBits,*_alphabet);

        // the preamble is set before the symbols to be transmitted
        if(_preamble.size()>0)
        {
            _symbols.resize(_preamble.rows(),_preamble.cols()+symbolsWithoutPreamble.cols());
            _symbols << _preamble,symbolsWithoutPreamble;
        }else
            _symbols = symbolsWithoutPreamble;

        // all the above symbols must be processed except those generated due to the smoothing
        _iLastSymbolVectorToBeDetected = _symbols.cols() - _nSmoothingSymbolsVectors;

        buildChannel();

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
		_noise = new SingleUserPowerProfileDependentNoise(_L,_channel->length(),*_powerProfile);

#ifdef EXPORT_REAL_DATA
            realSymbols = &_symbols;
            realChannel = _channel;
            realNoise = _noise;
#endif

        for(_iSNR=0;_iSNR<_SNRs.size();_iSNR++)
        {
            cout << COLOR_FRAME_NUMBER_SNR << "SNR = " << COLOR_NORMAL << _SNRs[_iSNR] << COLOR_FRAME_NUMBER_SNR << " [Frame " << COLOR_NORMAL << _iFrame << COLOR_FRAME_NUMBER_SNR << "]..." << COLOR_NORMAL << endl;

            // noise SNR is set
            _noise->setSNR(_SNRs[_iSNR],_alphabet->variance());

#ifdef PRINT_NOISE
			cout << "noise is" << endl;
			noise->print();
			cout << endl;
#endif

            // transmission
            _observations = _channel->transmit(_symbols,*_noise);

             addAlgorithms();

            // here the number of algoriths is known. So, the first iteration:
            if(_iFrame==0 && _iSNR==0)
                onlyOnce();

            // algorithms are executed
            for(_iAlgorithm=0;_iAlgorithm<_algorithms.size();_iAlgorithm++)
            {
                // if there is training sequence
                if(_trainSeqLength!=0)
                    _algorithms[_iAlgorithm]->run(_observations,_noise->variances(),_symbols.block(0,_preambleLength,_N,_trainSeqLength));
                // if there is NOT training sequence
                else
                    _algorithms[_iAlgorithm]->run(_observations,_noise->variances());

                _detectedSymbols = _algorithms[_iAlgorithm]->getDetectedSymbolVectors();

                _pe = computeSER(_symbols.block(0,_preambleLength,_N,_frameLength),_detectedSymbols,_isSymbolAccountedForDetection,_iBestPermutation,_bestPermutationSigns);

                beforeEndingAlgorithm();

                delete _algorithms[_iAlgorithm];
            }

            _algorithms.clear();
			
#ifdef STOP_AFTER_EACH_SNR
		getchar();
#endif
        } // for(int iSNR=0;iSNR<SNRs.size();iSNR++)

// only if the results are to be saved after every processed frame, we initialize the file pointer with a valid filename at each frame
// FIXME: the program is still trying to save the data all the time (the calls to write in the file are made anyway)
#ifdef SAVE_ALL_DATA_AFTER_PROCESSING_EACH_FRAME
		_f.open(_tempOutputFileName,ofstream::trunc);
// otherwise the file pointer only gets initialized for the end frame
#else
	if(_iFrame==_nFrames-1)
		_f.open(_tempOutputFileName,ofstream::trunc);
#endif

        beforeEndingFrame();
        
#ifdef SAVE_ALL_DATA_AFTER_PROCESSING_EACH_FRAME
		_f.close();
#else
	if(_iFrame==_nFrames-1)
		_f.close();
#endif

		// the temporal file is renamed as the final
		std::string mvCommand = string(MV_COMMAND) + string(" ") + string(_tempOutputFileName) + string(" ") + string(_outputFileName);
		int systemCommandReturn = system(mvCommand.c_str());
		cout << "moving operation returned " << systemCommandReturn << endl;

        // ---------------------------------------------------------

        _iFrame++;

        delete _channel;
		_channel = NULL;
		
        delete _noise;
		_noise = NULL;
		
#ifdef STOP_AFTER_EACH_FRAME
		getchar();
#endif
    } // while((iFrame<nFrames) && (!done))

    _overallPeMatrix *= 1.0/_iFrame;
    _overallMseMatrix *= 1.0/_iFrame;

    cout << "Overall SER:" << endl;
    Util::print(_overallPeMatrix);

    cout << "Overall MSE:" << endl;
    Util::print(_overallMseMatrix);
    
    // data saving
    _xmlFile << "</com>" << endl;
    _xmlFile.close();
}

void BaseSystem::onlyOnce()
{
    _overallPeMatrix = MatrixXd::Zero(_SNRs.size(),_algorithms.size());
    _presentFramePe = MatrixXd::Zero(_SNRs.size(),_algorithms.size());

    _overallMseMatrix = MatrixXd::Zero(_SNRs.size(),_algorithms.size());
    _presentFrameMSE = MatrixXd::Zero(_SNRs.size(),_algorithms.size());

    // Pe evolution
    for(uint i=0;i<_SNRs.size();i++)
    {
        _overallPeTimeEvolution[i] = MatrixXd(_algorithms.size(),_frameLength);
        _overallErrorsNumberTimeEvolution[i] = MatrixXi::Zero(_algorithms.size(),_frameLength);
    }

    // we fill the vector with the names of the algorithms
    for(uint iAlgorithm=0;iAlgorithm<_algorithms.size();iAlgorithm++)
        _algorithmsNames.push_back(_algorithms[iAlgorithm]->getName());

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    for(uint i=0;i<SNRs.size();i++)
        presentFrameMSEtimeEvolution[i] = MatrixXd(algorithms.size(),frameLength);
#endif

#ifdef KEEP_ALL_CHANNEL_ESTIMATIONS
  // channel estimations (iChannelMatrixRow,iChannelMatrixCol,iTimeInstant,iAlgorithm,iSNR)
  presentFrameChannelMatrixEstimations = std::vector<std::vector<std::vector<MatrixXd> > >(_SNRs.size(),std::vector<std::vector<MatrixXd> >(_algorithms.size()));
#endif
}

void BaseSystem::beforeEndingAlgorithm()
{
//     mse = algorithms[iAlgorithm]->MSE(channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1),permutations[_iBestPermutation],_bestPermutationSigns);
//     mse = algorithms[iAlgorithm]->MSE(channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1));
	
	// the channel matrices estimated by the algorithm are stored
	std::vector<MatrixXd> estimatedChannelMatrices = _algorithms[_iAlgorithm]->getEstimatedChannelMatrices();
	
	if(estimatedChannelMatrices.size()!=0)
	{
	  std::vector<MatrixXd> toCheckEstimatedChannelMatrices(estimatedChannelMatrices.begin()+_MSEwindowStart,estimatedChannelMatrices.end());
// 	  cout << "estimatedChannelMatrices.size() = " << estimatedChannelMatrices.size() << " toCheckEstimatedChannelMatrices.size() = " << toCheckEstimatedChannelMatrices.size() << endl;
// 	  cout << "channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1).size() = " << channel->range(preambleLength+MSEwindowStart,iLastSymbolVectorToBeDetected-1).size() << endl;

	  _mse = computeMSE(_channel->range(_preambleLength+_MSEwindowStart,_iLastSymbolVectorToBeDetected-1),toCheckEstimatedChannelMatrices);
	}else
	{
	  // it should be zero
	  _mse = computeMSE(_channel->range(_preambleLength+_MSEwindowStart,_iLastSymbolVectorToBeDetected-1),estimatedChannelMatrices);
	}
	
#ifdef MSE_TIME_EVOLUTION_COMPUTING
    VectorXd mseAlongTime = TransmissionUtil::MSEalongTime(algorithms[iAlgorithm]->getEstimatedChannelMatrices(),0,frameLength-1,channel->range(preambleLength,preambleLength+frameLength-1),0,frameLength-1);
    for(int ik=0;ik<frameLength;ik++)
        presentFrameMSEtimeEvolution[iSNR](iAlgorithm,ik) = mseAlongTime(ik);
#endif

#ifdef KEEP_ALL_CHANNEL_ESTIMATIONS
  // we get the channel matrices estimated by this algorithm
  vector<MatrixXd> thisAlgorithmEstimatedChannelMatrices = _algorithms[_iAlgorithm]->getEstimatedChannelMatrices();

  // if none, that meaning the algorithm does not performa channel matrix estimation,...
  if(thisAlgorithmEstimatedChannelMatrices.size()==0)
	// we generate a sequence of matrices
	thisAlgorithmEstimatedChannelMatrices = vector<MatrixXd>(_iLastSymbolVectorToBeDetected-_preambleLength,MatrixXd::Zero(_channel->channelCoefficientsMatrixRows(),_channel->channelCoefficientsMatrixCols()));
  else
	if(thisAlgorithmEstimatedChannelMatrices.size()!=static_cast<uint>(_iLastSymbolVectorToBeDetected-_preambleLength))
	  throw RuntimeException("BaseSystem::BeforeEndingAlgorithm: the number of channel matrices estimated by the algorithm is not the expected.");

  presentFrameChannelMatrixEstimations[_iSNR][_iAlgorithm] = thisAlgorithmEstimatedChannelMatrices;
#endif

    cout << COLOR_GREEN << _algorithms[_iAlgorithm]->getName() << COLOR_NORMAL << ": Pe = " << _pe << " , MSE = " << _mse << endl;

    // the error probability is accumulated
    _overallPeMatrix(_iSNR,_iAlgorithm) += _pe;
    _presentFramePe(_iSNR,_iAlgorithm) = _pe;

    // and the MSE
    _overallMseMatrix(_iSNR,_iAlgorithm) += _mse;
    _presentFrameMSE(_iSNR,_iAlgorithm) = _mse;

    // Pe evolution
    MatrixXd transmittedSymbols = _symbols.block(0,_preambleLength,_N,_frameLength);

    if(_detectedSymbols.rows()!=0)
    {
        for(int k=0;k<_frameLength;k++)
            for(uint iUser=0;iUser<_N;iUser++)
                if(_detectedSymbols(iUser,k)!=transmittedSymbols(iUser,k))
                    _overallErrorsNumberTimeEvolution[_iSNR](_iAlgorithm,k)++;
    }
}

void BaseSystem::beforeEndingFrame()
{
    // pe
    _peMatrices.push_back(_presentFramePe);
    Util::matricesVectorToOctaveFileStream(_peMatrices,"pe",_f);

    // MSE
    _MSEMatrices.push_back(_presentFrameMSE);
    Util::matricesVectorToOctaveFileStream(_MSEMatrices,"mse",_f);

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    MSEtimeEvolution.push_back(presentFrameMSEtimeEvolution);
    Util::matricesVectorsVectorToOctaveFileStream(MSEtimeEvolution,"MSEtimeEvolution",f);
#endif

#ifdef KEEP_ALL_CHANNEL_MATRICES
	_channelMatrices.push_back(_channel->range(_preambleLength,_iLastSymbolVectorToBeDetected-1));
#endif

//     for(uint iSNR=0;iSNR<SNRs.size();iSNR++)
//         for(uint i=0;i<algorithmsNames.size();i++)
//             for(int j=0;j<frameLength;j++)
//                 overallPeTimeEvolution[iSNR](i,j) = (double) overallErrorsNumberTimeEvolution[iSNR](i,j) / (double) (N*(iFrame+1));
//     Util::matricesVectorToOctaveFileStream(overallPeTimeEvolution,"peTimeEvolution",f);

    Util::scalarToOctaveFileStream(_iFrame+1,"nFrames",_f);

    Util::stringsVectorToOctaveFileStream(_algorithmsNames,"algorithmsNames",_f);
    Util::scalarToOctaveFileStream(_L,"L",_f);
    Util::scalarToOctaveFileStream(_N,"N",_f);
    Util::scalarToOctaveFileStream(_m,"m",_f);
    Util::scalarToOctaveFileStream(_frameLength,"frameLength",_f);
    Util::scalarToOctaveFileStream(_trainSeqLength,"trainSeqLength",_f);
    Util::scalarToOctaveFileStream(_d,"d",_f);
    Util::scalarToOctaveFileStream(_symbolsDetectionWindowStart,"symbolsDetectionWindowStart",_f);
    Util::scalarToOctaveFileStream(_MSEwindowStart,"MSEwindowStart",_f);
    Util::scalarsVectorToOctaveFileStream(_SNRs,"SNRs",_f);
    Util::matrixToOctaveFileStream(_preamble,"preamble",_f);
    Util::scalarToOctaveFileStream(_nSmoothingSymbolsVectors,"nSmoothingSymbolsVectors",_f);    
    Util::scalarToOctaveFileStream(_preambleLength,"preambleLength",_f);
    Util::scalarsVectorToOctaveFileStream(_mainSeeds,"mainSeeds",_f);
    Util::scalarsVectorToOctaveFileStream(_statUtilSeeds,"statUtilSeeds",_f);
	
#ifdef KEEP_ALL_CHANNEL_MATRICES
	Util::matricesVectorsVectorToOctaveFileStream(_channelMatrices,"channels",_f);
#else
	// only last channel is saved
    Util::matricesVectorToOctaveFileStream(_channel->range(_preambleLength,_iLastSymbolVectorToBeDetected),"channel",_f);
#endif

    Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*_channel).name())),"channelClass",_f);
    Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*_noise).name())),"noiseClass",_f);
    Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*this).name())),"systemClass",_f);

    if(_powerProfile!=NULL)
    {
        Util::scalarsVectorToOctaveFileStream(_powerProfile->tapsPowers(),"powerProfileVariances",_f);
        Util::stringsVectorToOctaveFileStream(vector<string>(1,string(typeid(*_powerProfile).name())),"powerProfileClass",_f);
    }
    
#ifdef KEEP_ALL_CHANNEL_ESTIMATIONS
	channelEstimations.push_back(presentFrameChannelMatrixEstimations);
	Util::matricesVectorsVectorsVectoresVectorToOctaveFileStream(channelEstimations,"channelEstimations",_f);
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

    if(sourceSymbols.rows()!= detectedSymbols.rows() || static_cast<uint>(detectedSymbols.rows())!= mask.size())
    {
        cout << "sourceSymbols.rows() = " << sourceSymbols.rows() << " detectedSymbols.rows() = " << detectedSymbols.rows() << " mask.size() = " << mask.size() << endl;
        throw RuntimeException("BaseSystem::computeSER: matrix row numbers differ.");
    }

    if(sourceSymbols.cols()!= detectedSymbols.cols() || static_cast<uint>(detectedSymbols.cols())!= mask[0].size())
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

    for(uint iPermut=0;iPermut<_permutations.size();iPermut++)
    {
        int permutationErrors = 0;
        
        for(uint iStream=0;iStream<_permutations[iPermut].size();iStream++)
        {
            iInput = _permutations[iPermut][iStream];
          
            int currentStreamErrorsInverting=0,currentStreamErrorsWithoutInverting=0;
            
            for(uint iTime=0;iTime<static_cast<uint> (sourceSymbols.cols());iTime++)
            {
                // if this symbol is not accounted for
                if(!mask[iStream][iTime])
                    continue;

                // if the symbols differ, an error happened...
                currentStreamErrorsWithoutInverting += (sourceSymbols(iStream,iTime) != detectedSymbols(iInput,iTime));
                
                // ...unless the symbol sign needs to be switched because of the ambiguity
                currentStreamErrorsInverting += (sourceSymbols(iStream,iTime) != _alphabet->opposite(detectedSymbols(iInput,iTime)));
                
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

            if(currentStreamErrorsWithoutInverting<currentStreamErrorsInverting)
            {
                permutationErrors += currentStreamErrorsWithoutInverting;
                thisPermutationSigns[iStream] = 1;
            }
            else
            {
                permutationErrors += currentStreamErrorsInverting;
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
    
	assert(nAccountedSymbols % _permutations.size() == 0);
    
    nAccountedSymbols /= _permutations.size();
    
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

  if(sourceSymbols.rows()!= detectedSymbols.rows() || static_cast<uint>(detectedSymbols.rows())!= mask.size())
  {
	  cout << "sourceSymbols.rows() = " << sourceSymbols.rows() << " detectedSymbols.rows() = " << detectedSymbols.rows() << " mask.size() = " << mask.size() << endl;
	  throw RuntimeException("BaseSystem::computeSERwithoutSolvingAmbiguity: matrix row numbers differ.");
  }

  if(sourceSymbols.cols()!= detectedSymbols.cols() || static_cast<uint>(detectedSymbols.cols())!= mask[0].size())
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

  for(uint iStream=0;iStream<_N;iStream++)
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
