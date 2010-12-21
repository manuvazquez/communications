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
#include "ISWCS10System.h"

#include <bashcolors.h>
#include <OneChannelOrderPerOutputSMCAlgorithm.h>

#include <PSPBasedSMCAlgorithm.h>

// #define USE_BESSEL_CHANNEL
#define USE_AR_CHANNEL

ISWCS10System::ISWCS10System()
 : ChannelOrderEstimationSystem()
{
// 	nSurvivors = 1;
	nSurvivors = 2;
    adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;

	_velocity = 50.0; // m/s
	_carrierFrequency = 2e9;
	_period = 1.0/500.0e3;


	// in order to use a Bessel channel (considering the Clarke model), the parameters of the AR process that the algorithms will consider are derived from those of the Clarke model
	double computedARprocessVariance;
	std::vector<double> computedARcoeffs = ARprocess::parametersFromYuleWalker(ARcoefficients.size(),_velocity,_carrierFrequency,_period,computedARprocessVariance);
	
// 	cout << "computedARprocessVariance = " << computedARprocessVariance << endl;
// 	Util::print(computedARcoeffs);
// 	cout << endl;
	

    _powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);

	if(adjustParticlesNumberFromSurvivors)
	{
		nParticles = (int)pow((double)_alphabet->length(),_N*(_m-1))*nSurvivors;
        cout << COLOR_WHITE << "Number of particles adjusted to " << COLOR_NORMAL << nParticles << endl;
    }

	if(adjustSurvivorsFromParticlesNumber)
	{
		cout << COLOR_WHITE << "Number of survivors adjusted from " << COLOR_NORMAL << nSurvivors;
		nSurvivors = int(ceil(double(nParticles)/pow(2.0,double(_N*(_m-1)))));
		cout << COLOR_WHITE << " to " << COLOR_NORMAL << nSurvivors << endl;
	}

	for(uint iChannelOrder=0;iChannelOrder<_candidateChannelOrders.size();iChannelOrder++)
	{
	  #if defined USE_AR_CHANNEL
		kalmanChannelEstimators.push_back(new KalmanEstimator(MatrixXd::Zero(1,_N*_candidateChannelOrders[iChannelOrder]),
															  MatrixXd::Ones(1,_N*_candidateChannelOrders[iChannelOrder]),_N,ARcoefficients,ARvariance));
	  #elif defined USE_BESSEL_CHANNEL
		  kalmanChannelEstimators.push_back(new KalmanEstimator(MatrixXd::Zero(1,_N*_candidateChannelOrders[iChannelOrder]),
															  MatrixXd::Ones(1,_N*_candidateChannelOrders[iChannelOrder]),_N,computedARcoeffs,computedARprocessVariance));
	  #else
		std::cout << "in ISWCS10System.cpp not defined a type channel to use" << std::endl;
		exit(1);
	  #endif
	}

    ResamplingCriterion resamplingCriterion(resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
	bestParticlesResamplingAlgorithm = new BestParticlesResamplingAlgorithm(resamplingCriterion);

	
	#if defined USE_AR_CHANNEL
	  _kalmanEstimatorForActualChannelOrder = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,ARcoefficients,ARvariance);
	  _kalmanEstimatorForMaximumChannelOrder = new KalmanEstimator(_channelOrderCoefficientsMeans[_iMaxChannelOrder],_channelOrderCoefficientsVariances[_iMaxChannelOrder],_N,ARcoefficients,ARvariance);
	#elif defined USE_BESSEL_CHANNEL
	  _kalmanEstimatorForActualChannelOrder = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,computedARcoeffs,computedARprocessVariance);
	  _kalmanEstimatorForMaximumChannelOrder = new KalmanEstimator(_channelOrderCoefficientsMeans[_iMaxChannelOrder],_channelOrderCoefficientsVariances[_iMaxChannelOrder],_N,computedARcoeffs,computedARprocessVariance);
	#else
	  std::cout << "in ISWCS10System.cpp not defined a type channel to use" << std::endl;
	  exit(1);
	#endif	
	
	
// 	// 1-3-1
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[1] = 3;

// 	// 3-1-3
// 	_subchannelOrders = std::vector<uint>(3,3);
// 	_subchannelOrders[1] = 1;

// 	// 3-3-3
// 	_subchannelOrders = std::vector<uint>(3,3);

// ---------------------------------

	// 4-1-1
	_subchannelOrders = std::vector<uint>(3,1);
	_subchannelOrders[0] = 4;

// 	// 4-3-2
// 	_subchannelOrders = std::vector<uint>(3);
// 	_subchannelOrders[0] = 4;_subchannelOrders[1] = 3;_subchannelOrders[2] = 2;

// 	// 4-4-4
// 	_subchannelOrders = std::vector<uint>(3,4);

// ----------------------------------

// 	// 2-1-1 (ambiguity problems!!)
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[0] = 2;

// 	// 2-2-1
// 	_subchannelOrders = std::vector<uint>(3,2);
// 	_subchannelOrders[2] = 1;

// 	// 4-3-2
// 	_subchannelOrders = std::vector<uint>(3,4);
// 	_subchannelOrders[2] = 3;
// 	_subchannelOrders[2] = 2;

// 	// 4-4-1
// 	_subchannelOrders = std::vector<uint>(3,4);
// 	_subchannelOrders[2] = 1;
}


ISWCS10System::~ISWCS10System()
{
	delete _powerProfile;

	for(uint iChannelOrder=0;iChannelOrder<_candidateChannelOrders.size();iChannelOrder++)
		delete kalmanChannelEstimators[iChannelOrder];

	delete withoutReplacementResamplingAlgorithm;
	delete bestParticlesResamplingAlgorithm;

	delete _kalmanEstimatorForActualChannelOrder;
	delete _kalmanEstimatorForMaximumChannelOrder;
}

void ISWCS10System::buildChannel()
{
//     channel = new ARchannel(N,L,m,symbols.cols(),ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
	_channel = new BesselChannel(_N,_L,_m,_symbols.cols(),_velocity,_carrierFrequency,_period,*_powerProfile);

	dynamic_cast<StillMemoryMIMOChannel*>(_channel)->setSubchannelOrders(_subchannelOrders);

// 	channel = new TimeInvariantChannel(N,L,m,symbols.cols(),powerProfile->generateChannelMatrix(randomGenerator));
}

void ISWCS10System::addAlgorithms()
{
	ChannelOrderEstimationSystem::addAlgorithms();

//  	_algorithms.push_back(new OneChannelOrderPerOutputSMCAlgorithm("OneChannelOrderPerOutputSMCAlgorithm",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,kalmanChannelEstimators,_preamble,_preamble.cols(),_d,nParticles,bestParticlesResamplingAlgorithm));

 	_algorithms.push_back(new PSPAlgorithm("PSPAlgorithm (known maximum suborder)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimatorForActualChannelOrder,_preamble,_d,_iLastSymbolVectorToBeDetected+_d,nSurvivors));

// // 	_algorithms.push_back(new PSPAlgorithm("PSPAlgorithm (maximum suborder within the set)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_candidateChannelOrders[_iMaxChannelOrder],_kalmanEstimatorForMaximumChannelOrder,_preamble,_candidateChannelOrders[_iMaxChannelOrder]-1,_iLastSymbolVectorToBeDetected+_candidateChannelOrders[_iMaxChannelOrder]-1,nSurvivors));

	_algorithms.push_back(new PSPBasedSMCAlgorithm("PSP based SMC algorithm (deterministic)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimatorForActualChannelOrder,_preamble,_d,nParticles,bestParticlesResamplingAlgorithm,_powerProfile->means(),_powerProfile->variances()));
	
//      _algorithms.push_back(new ViterbiAlgorithm("Viterbi (known channel)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (_channel)),_preamble,_d));
}

void ISWCS10System::beforeEndingFrame()
{
    ChannelOrderEstimationSystem::beforeEndingFrame();
    Util::scalarToOctaveFileStream(nSurvivors,"nSurvivors",_f);
	Util::scalarToOctaveFileStream(_velocity,"velocity",_f);
	Util::scalarsVectorToOctaveFileStream(_subchannelOrders,"subchannelOrders",_f);
}

// double ISWCS10System::computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns)
// {
// //   if(_symbolsDetectionWindowStart!=0)
// // 	throw RuntimeException("ISWCS10System::computeSER: this is only implemented when the starting point for SER computing is zero (the beginning of the frame).");
//   
//   if(detectedSymbols.rows()==0)
//   {
// 	_piecesBestPermuationIndexes = std::vector<uint>(_signChanges.size()-1,0);
// 	_piecesBestPermutationSigns = std::vector<std::vector<int> >(_signChanges.size()-1,std::vector<int>(_N,1));
// 	return -1.0;
//   }
// 
//   _piecesBestPermuationIndexes.clear();
//   _piecesBestPermutationSigns.clear();  
// 
//   double res = 0.0;
//   
// //   _signChanges = _channel->getInputsZeroCrossings(_preambleLength+_symbolsDetectionWindowStart,_frameLength-_preambleLength-_symbolsDetectionWindowStart);
//   _signChanges = _channel->getInputsZeroCrossings(_preambleLength+_symbolsDetectionWindowStart,_frameLength-_symbolsDetectionWindowStart);
//   //								^
//   //								|
//   // even though ultimately only symbol vectors from "symbolsDetectionWindowStart" will be taken into account for detection, we don't
//   // have to worry about it here, since the mask will take care of that. That being so, "_signChanges" actually contains all the sign changes
//   // that occur within the frame
// //   _signChanges = _channel->getInputsZeroCrossings(_preambleLength,_frameLength);
// 
//   uint overallNumberAccountedSymbols = 0;
//   
//   for(uint iSignChange=1;iSignChange<_signChanges.size();iSignChange++)
//   {
// 	  // we need to find out how many symbols are gonna be taken into account within this subframe
// 	  uint thisSubframeNumberAccountedSymbols = 0;
// 	  for(int i=0;i<_N;i++)
// 		for(uint j=_signChanges[iSignChange-1];j<_signChanges[iSignChange];j++)
// 		  thisSubframeNumberAccountedSymbols += mask[i][j-_preambleLength];
// 
// 	  overallNumberAccountedSymbols += thisSubframeNumberAccountedSymbols;
// 
// // 	  cout << "_frameLength = " << _frameLength << endl;
// // 	  cout << " _signChanges[iSignChange-1] = " << _signChanges[iSignChange-1] << endl;
// // 	  cout << " _signChanges[iSignChange] = " << _signChanges[iSignChange] << endl;
// // 	  cout << "_signChanges[iSignChange]-_signChanges[iSignChange-1] = " << _signChanges[iSignChange]-_signChanges[iSignChange-1] << endl;
// // 	  cout << " mask.rows() = " << mask.size() << " mask.cols() = " << mask[0].size() << endl;
// // 	  cout << "sourceSymbols.rows() = " << sourceSymbols.rows() << " sourceSymbols.cols() = " << sourceSymbols.cols() << endl;
// 	  
// 	  res += thisSubframeNumberAccountedSymbols*
// 				BaseSystem::computeSER(sourceSymbols.block(0,_signChanges[iSignChange-1]-_preambleLength,_N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
// 				detectedSymbols.block(0,_signChanges[iSignChange-1]-_preambleLength,_N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
// 				Util::block(mask,0,_signChanges[iSignChange-1]-_preambleLength,_N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
// 				iBestPermutation,bestPermutationSigns);
// 				
// 	  // we need to store which the best permutations were along with their corresponding signs
// 	  _piecesBestPermuationIndexes.push_back(iBestPermutation);
// 	  _piecesBestPermutationSigns.push_back(bestPermutationSigns);
//   }
// 
//   res /= overallNumberAccountedSymbols;
// 
//   return res;
// }