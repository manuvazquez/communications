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
#include <fstream>
#include <string>
#include <unistd.h>
#include <ctime>
#include <vector>

#include <defines.h>
#include <Alphabet.h>
#include <Bits.h>
#include <Modulator.h>
#include <Demodulator.h>
#include <Util.h>
#include <Octave.h>
#include <StatUtil.h>
#include <Random.h>
#include <MIMOChannel.h>
#include <ARchannel.h>
#include <Noise.h>
#include <ChannelDependentNoise.h>
#include <PowerProfileDependentNoise.h>
#include <NullNoise.h>
#include <Algorithm.h>

#include <bashcolors.h>

#include <rapidxml.hpp>

#define HOSTNAME_LENGTH 50
#define MV_COMMAND "mv"
#define LN_COMMAND "ln"

#define SYMBOLIC_LINK_NAME "last_res"

// to save all the generated channels (each one being a collection of channel matrices)
// #define KEEP_ALL_CHANNEL_MATRICES

// #define SAVE_ALL_SEEDS

// #define SAVE_NOISE_VARIANCES

extern bool __done,__saveData;

using namespace rapidxml;

class BaseSystem{
protected:
	bool _randomSeeds;
	bool _loadSeeds;
	bool _loadPerAlgorithmAndSNRseeds;
	
	uint32_t _mainRandomSeed;
	bool _mainRandomHaveStoredSample;
	double _mainRandomStoredSample;
	
	uint32_t _statUtilSeedToBeLoaded;
	
	uint32_t _perAlgorithmAndSNRstatUtilSeedToBeLoaded;
	
	bool _keepAllChannelEstimates,_keepAllDetectedSymbols;
	
#ifdef SAVE_CHANNEL_ESTIMATES_VARIANCES
	bool _saveChannelEstimatesVariances;
#endif
	
	std::string _noiseClassToBeInstantiated;
	std::string _channelClassToBeInstantiated;
	
    double _pe,_mse;
    uint _iSNR;
	uint _iFrame,_iAlgorithm;
	
	std::string _tmpResultsFile,_resultsFile;

    Alphabet *_alphabet;

    Noise *_noise;
    MatrixXd _observations;
	
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

	/**
	 * @brief SNRs to be processed
	 **/
	
    std::vector<int> _SNRs;

	/**
	 * @brief when SER computing starts (with respect to the beginning of the frame length)
	 **/
	uint _symbolsDetectionWindowStart;
	
	// when MSE computing starts (with respect to the beginning of the frame length)
	uint _MSEwindowStart;

    // a vector that will contain the names of the algorithms
    std::vector<std::string> _algorithmsNames;

    MatrixXd _preamble;

    // algorithms performing smoothing require symbol vector x_{frameLength:frameLength+d} in order to detect the last symbol vector
    uint _nSmoothingSymbolsVectors;

	/**
	 * @brief Bessel channel parameters
	 **/
	double _velocity,_carrierFrequency,_symbolRate,_T;

	/**
	 * @brief coefficients for the Auto-Regressive channel
	 **/
	std::vector<double> _ARcoefficients;
	
	/**
	 * @brief variance for the Auto-Regressive channel
	 **/
    double _ARvariance;
    
	// indicates wether or not a symbol must be taken into account for detection. It has a bool for every symbol WITHIN THE FRAME, and hence it doesn't include preamble symbols or smoothing symbols
    vector<vector<bool> > _isSymbolAccountedForDetection;
	
	// indicates wether or not a channel estimate must be taken into account for MSE computation. It has a bool for every channel estimate WITHIN THE FRAME, and hence it doesn't include channel estimates corresponding to the preamble or smoothing
	vector<bool> _isChannelEstimateAccountedForMSE;
    
    std::vector<std::vector<uint> > _permutations;
	uint _iBestPermutation;
	std::vector<int> _bestPermutationSigns;

    // matrices for results
    std::vector<MatrixXd> _peMatrices, _MSEMatrices;

    // matrices for accumulating the probabiliy of error (MSE) for all SNR's and all algorithms so that they can be printed when the program finishes (they are not saved)
    MatrixXd _overallPeMatrix,_overallMseMatrix;
    
    // matrices for accumulating the probabiliy of error (MSE) for all SNR's and all algorithms in order to save them
    MatrixXd _presentFramePe,_presentFrameMSE;

    // BER time evolution
    std::vector<MatrixXd> _overallPeTimeEvolution;
    std::vector<MatrixXi> _overallErrorsNumberTimeEvolution;

    // seeds
	std::vector<Random> _mainRandoms,_statUtilRandoms;
	
#ifdef SAVE_NOISE_VARIANCES
	std::vector<double> _presentFrameNoiseVariances;
	std::vector<std::vector<double> > _noiseVariances;
#endif

#ifdef SAVE_ALL_SEEDS
	std::vector<std::vector<std::vector<Random> > > _perAlgorithmAndSNRstatUtilRandoms;
	std::vector<std::vector<Random> > _thisFramePerAlgorithmAndSNRstatUtilRandoms;
#endif

std::vector<std::vector<std::vector<std::vector<MatrixXd> > > > _channelEstimates;
std::vector<std::vector<std::vector<MatrixXd> > >  _presentFrameChannelMatrixEstimates;

std::vector<std::vector<std::vector<MatrixXd> > > _allDetectedSymbols;
std::vector<std::vector<MatrixXd> > _presentFrameAllDetectedSymbols;

std::vector<std::vector<std::vector<std::vector<MatrixXd> > > > _channelEstimatesVariances;
std::vector<std::vector<std::vector<MatrixXd> > >  _presentFrameChannelEstimatesVariances;

#ifdef KEEP_ALL_CHANNEL_MATRICES
	std::vector<std::vector<MatrixXd> > _channelMatrices;
#endif

    Random _randomGenerator;

    MIMOChannel *_channel;

    std::vector<Algorithm *> _algorithms;

    MatrixXd _symbols;
    MatrixXd _detectedSymbols;

    ofstream _f;
    std::ifstream _parametersFile;
	
	char *_parameters;

    DelayPowerProfile *_powerProfile;
	
	bool _saveAtEveryFrame;

    virtual void addAlgorithms() = 0;
	virtual void buildSystemSpecificVariables() {}
	virtual void storeFrameResults();
	virtual void saveFrameResults();
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
	
	double computeSymbolVectorErrorRate(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask) const;
	
	//! It computes de Mean Square Error
	/*!
	  \param realChannelMatrices the sequence of actual channel matrices used during the transmission
	  \param detectedChannelMatrices the detected channel matrices
	  \return the computed MSE
	*/
	virtual double computeMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices,const std::vector<bool> &mask) const;
	
	virtual double computeMSE(const vector<MatrixXd> &realchannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices,const std::vector<bool> &mask,const vector<uint> &bestPermutation,const vector<int> &bestPermutationSigns) const;
	
	xml_node<>* get_child(xml_node<> *inputNode, std::string sNodeFilter);
	
	template<class T> void readParameterFromXML(xml_node<> *parentNode,std::string xmlName,T &parameter);
	template<class T> void readMultiValuedParameterFromXML(xml_node<> *parentNode,std::string xmlName,std::vector<T> &vector);
	
	xml_document<> _doc;
	
	/**
	 * @brief it builds and returns a pointer to a Noise object according to the content of the variable _noiseClassToBeInstantiated read from XML
	 *
	 * @return Noise*
	 **/
	virtual Noise *buildNoise() const;
	
	/**
	 * @brief it builds and returns a pointer to a MIMOChannel object according to the content of the variable  _channelClassToBeInstantiated read from XML (it might need to modify some class variables, e.g., in CDMA the activity)
	 *
	 * @return MIMOChannel*
	 **/
	virtual MIMOChannel *buildChannel();

	template<typename Func> double computeXXXwithoutSolvingAmbiguity(const MatrixXd &sourceSymbols,const MatrixXd &obtainedSymbols,const vector<vector<bool> > &mask, Func func) const
	{
		assert( (sourceSymbols.rows() == obtainedSymbols.rows()) && (static_cast<uint>(obtainedSymbols.rows())== mask.size()) );
		assert( (sourceSymbols.cols()== obtainedSymbols.cols()) && (static_cast<uint>(obtainedSymbols.cols())== mask[0].size()) );

		uint nSymbolsRows = obtainedSymbols.rows();

		double error = 0.0;

		uint nAccountedSymbols = 0;

		for(uint iStream=0;iStream<nSymbolsRows;iStream++)
		{
			for(uint iTime=0;iTime<sourceSymbols.cols();iTime++)
			{
				// if this symbol is not accounted for
				if(!mask[iStream][iTime])
					continue;

				// if the symbols differ, an error happened...
				error += func(sourceSymbols(iStream,iTime),obtainedSymbols(iStream,iTime));
				
				nAccountedSymbols++;
			}              
		} // for(uint iStream=0;iStream<permutations[iPermut].size();iStream++)

		// if all the symbols were masked
		if(nAccountedSymbols==0)
			return 0.0;
		else
			return error/(double)(nAccountedSymbols);
	}
	
	double computeSymbolsMSEwithoutSolvingAmbiguity(const MatrixXd &sourceSymbols,const MatrixXd &estimatedSymbols,const vector<vector<bool> > &mask) const
	{
		return computeXXXwithoutSolvingAmbiguity(sourceSymbols,estimatedSymbols,mask,[] (const double a, const double b) {return (a-b)*(a-b);});
	}
	
	double computeSERwithoutSolvingAmbiguity(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask) const
	{
		return computeXXXwithoutSolvingAmbiguity(sourceSymbols,detectedSymbols,mask,[] (const double a, const double b) {return (a==b)?0.0:1.0; });
	}
	
public:
    BaseSystem();
    virtual ~BaseSystem();

    void simulate();
};

#endif
