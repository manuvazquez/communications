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

#include <lapackpp/gmd.h>
#include <lapackpp/gmi.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>
#include <lapackpp/sybmd.h>
#include <lapackpp/sybfd.h>

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
#include <Algorithm.h>
#include <TransmissionUtil.h>

#define HOSTNAME_LENGTH 50

#define MSE_TIME_EVOLUTION_COMPUTING

#ifdef EXPORT_REAL_DATA
    MIMOChannel *realChannel;
    tMatrix *realSymbols;
    Noise *realNoise;
#endif

extern bool __done;

class BaseSystem{
protected:
    double pe,mse;
    uint iSNR;
    int lastSymbolVectorInstant;

    // GLOBAL PARAMETERS
    int nFrames,L,N,K,m,d,trainSeqLength,preambleLength;
    char outputFileName[HOSTNAME_LENGTH+4];

    Alphabet *alphabet;

	Noise *ruido;
	tMatrix observaciones;

    // SNRs to be processed
    std::vector<int> SNRs;

    // BER and MSE computing
    int BERwindowStart,MSEwindowStart;

    // a vector that will contain the names of the algorithms
    std::vector<std::string> algorithmsNames;

    tMatrix preamble;

    int nSmoothingSymbolsVectors,nSmoothingBitsVectors;

    std::vector<std::vector<uint> > permutations;

    // matrices for results
    vector<tMatrix> peMatrices, MSEMatrices;

    tMatrix overallPeMatrix,overallMseMatrix,presentFramePe,presentFrameMSE;

    // BER time evolution
    std::vector<tMatrix> overallPeTimeEvolution;
    std::vector<LaGenMatInt> overallErrorsNumberTimeEvolution;

    // seeds
    std::vector<uint32_t> mainSeeds,statUtilSeeds;
    std::vector<LaGenMatLongInt> beforeRunStatUtilSeeds;
    LaGenMatLongInt presentFrameStatUtilSeeds;

#ifdef MSE_TIME_EVOLUTION_COMPUTING
    vector<tMatrix> presentFrameMSEtimeEvolution;
    vector<vector<tMatrix> > MSEtimeEvolution;
#endif

    Random randomGenerator;

    MIMOChannel *channel;

    std::vector<Algorithm *> algorithms;

    tMatrix symbols;
    tMatrix detectedSymbols;

    ofstream f;

    DelayPowerProfile *powerProfile;

    virtual void AddAlgorithms() = 0;
    virtual void BuildChannel() = 0;
    virtual void BeforeEndingFrame(int iFrame);
    virtual void BeforeEndingAlgorithm(int iAlgorithm);
    /**
     * Actions performed only once (first SNR, firs frame)
     */
    virtual void OnlyOnce();
public:
    BaseSystem();

    virtual ~BaseSystem();

    void Simulate();

};

#endif
