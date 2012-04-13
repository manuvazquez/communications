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
#include <TimeInvariantChannel.h>

// FIXME: do something to avoid the defines below: it's such a pain in the ass having to remember to change this when switching to a different kind of channel in the code

// #define USE_AR_CHANNEL
#define USE_BESSEL_CHANNEL
// #define USE_TIME_INVARIANT_CHANNEL

ISWCS10System::ISWCS10System()
 : ChannelOrderEstimationSystem()
{
	nSurvivors = 2;

	_nParticles = 128;

    adjustSurvivorsFromParticlesNumber = true;
    adjustParticlesNumberFromSurvivors = false;

	// in order to use a Bessel channel (considering the Clarke model), the parameters of the AR process the algorithms will consider
	// are derived from those of the Clarke model
	double computedARprocessVariance;
	std::vector<double> computedARcoeffs = ARprocess::parametersFromYuleWalker(_ARcoefficients.size(),_velocity,_carrierFrequency,_period,computedARprocessVariance);
	
	// by default, AR channel is assumed
	std::vector<double> kalmanEstimatorARcoeffs = _ARcoefficients;
	double kalmanEstimatorVariance =  _ARvariance;
	
	#if defined USE_AR_CHANNEL
		// it's ok: by default, above it is assumed that the channel is AR
		
		std::cout << COLOR_MAROON << "assuming AR channel in generating Kalman estimators" << COLOR_NORMAL << std::endl;
	#elif defined USE_BESSEL_CHANNEL
// 		kalmanEstimatorARcoeffs = computedARcoeffs;
// 		kalmanEstimatorVariance = computedARprocessVariance;
// 		std::cout << COLOR_MAROON;
// 		std::cout << "AR process parameters computed using Yule-Walker:" << std::endl;
// 		std::cout << "\t AR coefficientes: ";
// 		Util::print(computedARcoeffs);
// 		std::cout << std::endl << "\t AR variance: " << computedARprocessVariance << std::endl;
// 		std::cout << COLOR_NORMAL;
	#elif defined USE_TIME_INVARIANT_CHANNEL
		kalmanEstimatorARcoeffs = std::vector<double>(ARcoefficients.size(),0);
		kalmanEstimatorARcoeffs[0] = 1.0;
		kalmanEstimatorVariance = 0;
		std::cout << COLOR_MAROON << "assuming time invariant channel in generating Kalman estimators" << COLOR_NORMAL << std::endl;
	#else
		std::cout << "ISWCS10System::ISWCS10System: channel type not #defined" << std::endl;
		exit(1);
	#endif

    _powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);
	
	adjustParticlesSurvivors(_nParticles,nSurvivors,adjustParticlesNumberFromSurvivors,adjustSurvivorsFromParticlesNumber);

	for(uint iChannelOrder=0;iChannelOrder<_candidateChannelOrders.size();iChannelOrder++)
	{
		kalmanChannelEstimators.push_back(new KalmanEstimator(
											MatrixXd::Zero(1,_N*_candidateChannelOrders[iChannelOrder]),
											MatrixXd::Ones(1,_N*_candidateChannelOrders[iChannelOrder]),
											_N,kalmanEstimatorARcoeffs,kalmanEstimatorVariance));
											
		kalmanWholeChannelEstimators.push_back(new KalmanEstimator(
											MatrixXd::Zero(_L,_N*_candidateChannelOrders[iChannelOrder]),
											MatrixXd::Ones(_L,_N*_candidateChannelOrders[iChannelOrder]),
											_N,kalmanEstimatorARcoeffs,kalmanEstimatorVariance));
	}

    ResamplingCriterion resamplingCriterion(_resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
	bestParticlesResamplingAlgorithm = new BestParticlesResamplingAlgorithm(resamplingCriterion);

	
	_kalmanEstimatorForActualChannelOrder = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,kalmanEstimatorARcoeffs,kalmanEstimatorVariance);
	_kalmanEstimatorForMaximumChannelOrder = new KalmanEstimator(
											_channelOrderCoefficientsMeans[_iMaxChannelOrder],
											_channelOrderCoefficientsVariances[_iMaxChannelOrder],
											_N,kalmanEstimatorARcoeffs,kalmanEstimatorVariance);
	
// 	// 1-3-1
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[1] = 3;

// 	// 3-1-3
// 	_subchannelOrders = std::vector<uint>(3,3);
// 	_subchannelOrders[1] = 1;

// 	// 3-3-3
// 	_subchannelOrders = std::vector<uint>(3,3);

// ---------------------------------

// 	// 4-1-1
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[0] = 4;

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

	// 4-4-1
	_subchannelOrders = std::vector<uint>(3,4);
	_subchannelOrders[2] = 1;
}


ISWCS10System::~ISWCS10System()
{
	delete _powerProfile;

	for(uint iChannelOrder=0;iChannelOrder<_candidateChannelOrders.size();iChannelOrder++)
	{
		delete kalmanChannelEstimators[iChannelOrder];
		delete kalmanWholeChannelEstimators[iChannelOrder];
	}

	delete withoutReplacementResamplingAlgorithm;
	delete bestParticlesResamplingAlgorithm;

	delete _kalmanEstimatorForActualChannelOrder;
	delete _kalmanEstimatorForMaximumChannelOrder;
}

void ISWCS10System::buildSystemSpecificVariables()
{
// 	channel = new ARchannel(N,L,m,symbols.cols(),ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
// 	_channel = new TimeInvariantChannel(_N,_L,_m,_symbols.cols(),_powerProfile->generateChannelMatrix(_randomGenerator));
	_channel = new BesselChannel(_N,_L,_m,_symbols.cols(),_velocity,_carrierFrequency,_period,*_powerProfile);

	dynamic_cast<StillMemoryMIMOChannel*>(_channel)->setSubchannelOrders(_subchannelOrders);
}

void ISWCS10System::addAlgorithms()
{
	ChannelOrderEstimationSystem::addAlgorithms();

 	_algorithms.push_back(new OneChannelOrderPerOutputSMCAlgorithm("MLSD-m",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,kalmanChannelEstimators,_preamble,_preamble.cols(),_d,_nParticles,bestParticlesResamplingAlgorithm));

//  	_algorithms.push_back(new PSPAlgorithm("PSP",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimatorForActualChannelOrder,_preamble,_d,_iLastSymbolVectorToBeDetected+_d,nSurvivors));

// 	// not used for the journal paper
// 	_algorithms.push_back(new PSPAlgorithm("PSPAlgorithm (maximum suborder among the receiving antennas)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_candidateChannelOrders[_iMaxChannelOrder],_kalmanEstimatorForMaximumChannelOrder,_preamble,_candidateChannelOrders[_iMaxChannelOrder]-1,_iLastSymbolVectorToBeDetected+_candidateChannelOrders[_iMaxChannelOrder]-1,nSurvivors));

	_algorithms.push_back(new PSPBasedSMCAlgorithm("G-PSP",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,_kalmanEstimatorForActualChannelOrder,_preamble,_d,_nParticles,bestParticlesResamplingAlgorithm,_powerProfile->means(),_powerProfile->variances()));
	
     _algorithms.push_back(new ViterbiAlgorithm("Viterbi (known channel)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (_channel)),_preamble,_d));
	 
	_algorithms.push_back(new MLSDmAlgorithm("MLSD-m (single channel order)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,kalmanWholeChannelEstimators,_preamble,_preamble.cols(),_d,_nParticles,bestParticlesResamplingAlgorithm,_ARcoefficients[0],_firstSampledChannelMatrixVariance,_ARvariance));
}

void ISWCS10System::saveFrameResults()
{
    ChannelOrderEstimationSystem::saveFrameResults();
    Octave::toOctaveFileStream(nSurvivors,"nSurvivors",_f);
	Octave::toOctaveFileStream(_subchannelOrders,"subchannelOrders",_f);
}
