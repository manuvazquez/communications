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
#ifndef BASESYSTEM_H
#define BASESYSTEM_H

/**
    @author Manu <manu@rustneversleeps>
*/


#include <iostream>
#include <string>
#include <unistd.h>
#include <ctime>
#include <vector>

#include <Alphabet.h>
#include <Bits.h>
#include <Modulator.h>
#include <Demodulator.h>
#include <Util.h>
#include <StatUtil.h>
#include <Random.h>
#include <MIMOChannel.h>
#include <Noise.h>
#include <ChannelDependentNoise.h>
#include <PowerProfileDependentNoise.h>
#include <NullNoise.h>
#include <Algorithm.h>
#include <TransmissionUtil.h>

#define HOSTNAME_LENGTH 50
#define MV_COMMAND "mv"
#define LN_COMMAND "ln"

#define SYMBOLIC_LINK_NAME "last_res"

// #define MSE_TIME_EVOLUTION_COMPUTING

// to save all the generated channels (each one being a collection of channel matrices)
#define KEEP_ALL_CHANNEL_MATRICES

#define KEEP_ALL_CHANNEL_ESTIMATIONS

#define SAVE_ALL_SEEDS

extern bool __done;

class BaseSystem{
protected:
    double _pe,_mse;
    uint _iSNR;
	uint _iFrame,_iAlgorithm;
	
	/**
	 * @brief the algorithms will stop the processing (and hence, the detection of the symbols transmitted) here (this is usually the last time instant for which enough observations are available to perform smoothing)
	 **/
    uint _iLastSymbolVectorToBeDetected;

    // GLOBAL PARAMETERS

	
	/**
	 * @brief number of outputs
	 **/
    uint _L;

	uint _frameLength,_nBitsGenerated,_m,_d,_trainSeqLength,_preambleLength;
	
	/**
	 * @brief number of frames to be simulated
	 **/
	uint _nFrames;
	
	/**
	 * @brief number of inputs
	 **/
	uint _N;
	
	std::string _tmpResultsFile,_resultsFile;

    Alphabet *_alphabet;

    Noise *_noise;
    MatrixXd _observations;

    // SNRs to be processed
    std::vector<int> _SNRs;

	/**
	 * @brief when SER computing starts (with respect to the beginning of the frame length)
	 **/
	uint _symbolsDetectionWindowStart;
	
	// when MSE computing starts (with respect to the beginning of the frame length)
	int _MSEwindowStart;

    // a vector that will contain the names of the algorithms
    std::vector<std::string> _algorithmsNames;

    MatrixXd _preamble;

    // algorithms performing smoothing require symbol vector x_{frameLength:frameLength+d} in order to detect the last symbol vector
    uint _nSmoothingSymbolsVectors;
    
	/*
	  indicates wether or not a symbol must be taken into account for detection. NOTE: this only has a bool for every information symbol
	  Hence,it doesn't include preamble symbols or smoothing symbols.
	*/
    vector<vector<bool> > _isSymbolAccountedForDetection;
    
    std::vector<std::vector<uint> > _permutations;
	uint _iBestPermutation;
	std::vector<int> _bestPermutationSigns;

    // matrices for results
    std::vector<MatrixXd> _peMatrices, _MSEMatrices;

    // matrices for accumulating the probabiliy of error (MSE) for all SNR's and all algorithms...
    // ...so that they can be printed when the program finishes (they are not saved)
    MatrixXd _overallPeMatrix,_overallMseMatrix;
    
    // matrices for accumulating the probabiliy of error (MSE) for all SNR's and all algorithms in order to save them
    MatrixXd _presentFramePe,_presentFrameMSE;

    // BER time evolution
    std::vector<MatrixXd> _overallPeTimeEvolution;
    std::vector<MatrixXi> _overallErrorsNumberTimeEvolution;

    // seeds
    std::vector<uint32_t> _mainSeeds,_statUtilSeeds;

#ifdef SAVE_ALL_SEEDS
	std::vector<std::vector<std::vector<uint32_t> > > /*_perAlgorithmAndSNRmainSeeds,*/_perAlgorithmAndSNRstatUtilSeeds;
	std::vector<std::vector<uint32_t> > /*_thisFramePerAlgorithmAndSNRmainSeeds,*/_thisFramePerAlgorithmAndSNRstatUtilSeeds;
#endif

#ifdef KEEP_ALL_CHANNEL_ESTIMATIONS
	std::vector<std::vector<std::vector<std::vector<MatrixXd> > > > _channelEstimations;
	std::vector<std::vector<std::vector<MatrixXd> > >  _presentFrameChannelMatrixEstimations;
#endif

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    vector<MatrixXd> presentFrameMSEtimeEvolution;
    vector<vector<MatrixXd> > MSEtimeEvolution;
#endif

#ifdef KEEP_ALL_CHANNEL_MATRICES
	std::vector<std::vector<MatrixXd> > _channelMatrices;
#endif

    Random _randomGenerator;

    MIMOChannel *_channel;

    std::vector<Algorithm *> _algorithms;

    MatrixXd _symbols;
    MatrixXd _detectedSymbols;

    ofstream _f,_xmlFile;

    DelayPowerProfile *_powerProfile;
	
	bool _saveAtEveryFrame;

    virtual void addAlgorithms() = 0;
    virtual void buildSystemSpecificVariables() = 0;
    virtual void beforeEndingFrame();
    virtual void beforeEndingAlgorithm();
	
    /*!
     * Actions performed only once (first SNR, first frame)
     */
    virtual void onlyOnce();

	//! It computes de Symbol Error Rate
	/*!
	  \param sourceSymbols the sequence of symbols actually transmitted
	  \param detectedSymbols the sequence of symbols detected
	  \param mask a matrix of bools indicating which symbols are to be taken into account in the calculus (in case not all of them are relevant)
	  \param iBestPermutation an int which will serve to store the best permutation
	  \param bestPermutationSigns a vector of unsigned ints where the signs of the best permutation will be stored
	  \return the computed probability
	*/
	virtual double computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns);
	
	double computeSERwithoutSolvingAmbiguity(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask) const;
	
// 	virtual double computeBER(const Bits &sourceBits,const Bits &detectedBits,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns);

	//! It computes de Mean Square Error
	/*!
	  \param realChannelMatrices the sequence of actual channel matrices used during the transmission
	  \param detectedChannelMatrices the detected channel matrices
	  \return the computed MSE
	*/
	virtual double computeMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices) const;
	
	virtual double computeMSE(const vector<MatrixXd> &realchannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices,const vector<uint> &bestPermutation,const vector<int> &bestPermutationSigns) const;

public:
    BaseSystem();
    virtual ~BaseSystem();

    void simulate();
};

#endif
