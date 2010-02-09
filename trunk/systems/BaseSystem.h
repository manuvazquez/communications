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
#include <sys/time.h>
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

// #define MSE_TIME_EVOLUTION_COMPUTING

#define KEEP_ALL_CHANNEL_MATRICES

extern bool __done;

class BaseSystem{
protected:
    double pe,mse;
    uint iSNR;
    int iLastSymbolVectorToBeDetected;

    // GLOBAL PARAMETERS
    int nFrames,L,N,frameLength,nBitsGenerated,m,d,trainSeqLength,preambleLength;
    char outputFileName[HOSTNAME_LENGTH+4];

    Alphabet *alphabet;

    Noise *noise;
    MatrixXd observations;

    // SNRs to be processed
    std::vector<int> SNRs;

    // BER and MSE computing
    int MSEwindowStart,symbolsDetectionWindowStart;

    // a vector that will contain the names of the algorithms
    std::vector<std::string> algorithmsNames;

    MatrixXd preamble;

    // algorithms performing smoothing require symbol vector x_{frameLength:frameLength+d} in order to detect the last symbol vector
    int nSmoothingSymbolsVectors;
    
    vector<vector<bool> > isSymbolAccountedForDetection;

    std::vector<std::vector<uint> > permutations;

    // matrices for results
    vector<MatrixXd> peMatrices, MSEMatrices;

    // matrices for accumulating the probabiliy of error (MSE) for all SNR's and all algorithms...
    // ...so that they can be printed when the program finishes (they are not saved)
    MatrixXd overallPeMatrix,overallMseMatrix;
    
    // matrices for accumulating the probabiliy of error (MSE) for all SNR's and all algorithms in order to save them
    MatrixXd presentFramePe,presentFrameMSE;

    // BER time evolution
    std::vector<MatrixXd> overallPeTimeEvolution;
    std::vector<MatrixXi> overallErrorsNumberTimeEvolution;

    // seeds
    std::vector<uint32_t> mainSeeds,statUtilSeeds;
//     std::vector<LaGenMatLongInt> beforeRunStatUtilSeeds;
//     LaGenMatLongInt presentFrameStatUtilSeeds;

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    vector<MatrixXd> presentFrameMSEtimeEvolution;
    vector<vector<MatrixXd> > MSEtimeEvolution;
#endif

#ifdef KEEP_ALL_CHANNEL_MATRICES
	std::vector<std::vector<MatrixXd> > channelMatrices;
#endif

    Random randomGenerator;

    MIMOChannel *channel;

    std::vector<Algorithm *> algorithms;

    MatrixXd symbols;
    MatrixXd detectedSymbols;

    ofstream f,xmlFile;

    DelayPowerProfile *powerProfile;

    virtual void AddAlgorithms() = 0;
    virtual void BuildChannel() = 0;
    virtual void BeforeEndingFrame(int iFrame);
    virtual void BeforeEndingAlgorithm(int iAlgorithm);
	
    /*!
     * Actions performed only once (first SNR, first frame)
     */
    virtual void OnlyOnce();

	//! It computes de symbol error rate
	/*!
	  \param sourceSymbols the sequence of symbols actually transmitted
	  \param detectedSymbols the sequence of symbols detected
	  \param mask a matrix of bools indicating which symbols are to be taken into account in the calculus (in case not all of them are relevant)
	  \param permutations a matrix containing the permutations within the symbol vectors that are to be explored
	  \return the computed probability
	*/
	virtual double computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask);

public:
    BaseSystem();

    virtual ~BaseSystem();

    void Simulate();

};

#endif
