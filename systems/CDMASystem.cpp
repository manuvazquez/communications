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
#include "CDMASystem.h"

#include <CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge.h>
#include <KnownFlatChannelOptimalAlgorithm.h>
#include <KnownFlatChannelAndActiveUsersOptimalAlgorithm.h>
#include <UnknownActiveUsersLinearFilterBasedSMCAlgorithm.h>
#include <OldUnknownActiveUsersLinearFilterBasedSMCAlgorithm.h>
#include <CDMAunknownActiveUsersSISopt.h>
#include <TimeInvariantChannel.h>
#include <MultiuserCDMAchannel.h>
#include <ViterbiAlgorithmWithAprioriProbabilities.h>
#include <PSPAlgorithmWithAprioriProbabilities.h>

#include <math.h>

#include <bashcolors.h>
#include <defines.h>
#include <SingleUserChannelDependentNoise.h>

#define PRINT_CODES_INFO
// #define PRINT_ACTIVITY_SAMPLING

// #define PRINT_INFO

// #define DEBUG_SER
// #define DEBUG_ACTIVITY_DETECTION_ERROR_RATE

// #define DEBUG

CDMASystem::CDMASystem(): SMCSystem()
,_userPersistenceProb(0.99),_newActiveUserProb(0.01),_userPriorProb(0.5)	// <- a reasonable model
// ,userPersistenceProb(0.8),newActiveUserProb(0.2),userPriorProb(1.0)
// ,_userPersistenceProb(1.0),_newActiveUserProb(0.2),_userPriorProb(1.0)		// all the users are active all the time
,_usersActivityPdfs(_N,UsersActivityDistribution(_userPersistenceProb,_newActiveUserProb,_userPriorProb))
,_iUserOfInterest(0u)
{
    if (_m!=1)
        throw RuntimeException("CDMASystem::CDMASystem: channel is not flat.");
	
 	_minSignalToInterferenceRatio = -30;
	//_minSignalToInterferenceRatio = -50;

	// first user starts transmitting something
	_usersActivityPdfs[0].setApriori(1.0);

	_nSurvivors = 2;

    // AR process parameters
    _ARcoefficients = vector<double>(2);
    _ARcoefficients[0] = 0.59999;
    _ARcoefficients[1] = 0.39999;
    _ARvariance=0.0001;    
    
//     // AR process parameters
//     ARcoefficients = vector<double>(1);
//     ARcoefficients[0] = 0.99999;
//     ARvariance=0.0001;
    
    // a flat power profile is generated. Notice:
    //      i) that m should be 1, otherwise an exception would have been thrown
    //     ii) we only need to generate a coefficient per user, i.e., a 1xN vector
// 	FlatPowerProfile::setCoefficientsMean(5.0);
    _powerProfile = new FlatPowerProfile(1,_N,_m,1.0);

	// bessel channel parameters
//     _velocity = 180/3.6; // (m/s)
    _velocity = 50/3.6; // (m/s)
    _carrierFrequency = 2e9; // (Hz)
    _symbolRate = 500e3; // (Hz)

    _T = 1.0/_symbolRate; // (s)
    

    // the parameters obtained from the Yule-Walker equations are used for the channel estimator
// 	_ARcoefficients = ARprocess::parametersFromYuleWalker(2,_velocity,_carrierFrequency,_T,_ARvariance);
    
	_cdmaKalmanEstimator = NULL;
	_cdmaKnownChannelChannelMatrixEstimator = NULL;

	_mmseDetector = new MMSEDetector(_L,_N,_alphabet->variance(),_N);
	
    // adjusting the number of particles from that of the survivors or the other way around
    _adjustSurvivorsFromParticlesNumber = false;
    _adjustParticlesNumberFromSurvivors = true;
	
    // check the adjustments for particle and survivor numbers
    if(_adjustParticlesNumberFromSurvivors && _adjustSurvivorsFromParticlesNumber)
        throw RuntimeException("CDMASystem::CDMASystem: \"adjustParticlesNumberFromSurvivors\" and \"adjustSurvivorsFromParticlesNumber\" shouldn't be true at the same time.");

    if(_adjustParticlesNumberFromSurvivors)
    {
	  cout << "Number of particles adjusted from " << nParticles;
	  
	  // the number of particles must be the number of states of the Viterbi/PSP algorithm times that of survivors
	  nParticles = (uint)pow((double)_alphabet->length()+1,_N)*_nSurvivors;
	  
	  cout << " to " << nParticles << endl;
    }

    if(_adjustSurvivorsFromParticlesNumber)
    {
	  cout << "Number of survivors adjusted from " << _nSurvivors;
	  _nSurvivors = uint(ceil(double(nParticles)/pow((double)_alphabet->length()+1,double(_N))));
	  cout << " to " << _nSurvivors << endl;
    }

	// in order to compute the BER/MSE/... of the user of interest only, a "dummy" permutation (a permutation for a set of 1 element) is generated
	_permutations = std::vector<std::vector<uint> >(1,std::vector<uint>(1,0));
	
	resetFramePieces();
		
	_everyFrameUsersActivity.reserve(_nFrames);
	_everyFrameSpreadingCodes.reserve(_nFrames);
	_peActivityDetectionFrames.reserve(_nFrames);
}


CDMASystem::~CDMASystem()
{
    delete _powerProfile;
    delete _cdmaKalmanEstimator;
	delete _cdmaKnownChannelChannelMatrixEstimator;
    delete _mmseDetector;
}

void CDMASystem::addAlgorithms()
{
	// ...the same for an estimator that knows the codes if these also change across frames
	delete _cdmaKalmanEstimator;
	_cdmaKalmanEstimator = new CDMAKalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_ARcoefficients,_ARvariance,_spreadingCodes);
	
    delete _cdmaKnownChannelChannelMatrixEstimator;
    _cdmaKnownChannelChannelMatrixEstimator = new CDMAKnownChannelChannelMatrixEstimator(_channel,_preambleLength,_N,_spreadingCodes);
	
    _algorithms.push_back(new KnownFlatChannelOptimalAlgorithm ("CDMA optimal with known channel BUT no knowledge of users activity probabilities",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,*_channel,_preambleLength));

    _algorithms.push_back(new KnownFlatChannelAndActiveUsersOptimalAlgorithm ("CDMA optimal (known channel and active users)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,*_channel,_preambleLength,_usersActivity));
     
    _algorithms.push_back(new CDMAunknownActiveUsersSISopt ("CDMA SIS-opt",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_preamble,_d,nParticles,algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_usersActivityPdfs));
	
    _algorithms.push_back(new CDMAunknownActiveUsersSISopt ("CDMA SIS-opt (known channel)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKnownChannelChannelMatrixEstimator,_preamble,_d,nParticles,algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_usersActivityPdfs));

//     _algorithms.push_back(new UnknownActiveUsersLinearFilterBasedSMCAlgorithm ("CDMA SIS Linear Filters",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_mmseDetector,_preamble,_d,nParticles,algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_usersActivityPdfs));
	
	_algorithms.push_back(new OldUnknownActiveUsersLinearFilterBasedSMCAlgorithm ("Old CDMA SIS Linear Filters",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_mmseDetector,_preamble,_d,nParticles,algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_usersActivityPdfs));

	_algorithms.push_back(new ViterbiAlgorithmWithAprioriProbabilities("Viterbi with a priori probabilities (known channel)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (_channel)),_preamble,_d,_usersActivityPdfs));
	
	_algorithms.push_back(new PSPAlgorithmWithAprioriProbabilities("PSP",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_preamble,_d,_iLastSymbolVectorToBeDetected+_d,_nSurvivors,_usersActivityPdfs));
	
	_algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_preamble,_symbols));
}

void CDMASystem::buildSystemSpecificVariables()
{
	// when users are not transmitting, their symbols are zero
	_usersActivity = vector<vector<bool> >(_symbols.rows(),vector<bool>(_frameLength));
    
    for(uint iUser=0;iUser<static_cast<uint>(_symbols.rows());iUser++)
    {
		// at the first time instant the prior probability is used to decide which users are active
		_usersActivity[iUser][_trainSeqLength] = _usersActivityPdfs[iUser].sampleFromPrior();
		
#ifdef PRINT_ACTIVITY_SAMPLING
		cout << "user " << iUser << ": " << _usersActivity[iUser][trainSeqLength] << endl;
#endif
		// when users are not transmitting, their symbols are zero
		_symbols(iUser,_preambleLength+_trainSeqLength) = double(_usersActivity[iUser][_trainSeqLength])*_symbols(iUser,_preambleLength+_trainSeqLength);
		
		// the symbol is accounted for detection only if the corresponding user is active
		_isSymbolAccountedForDetection[iUser][_trainSeqLength] = _usersActivity[iUser][_trainSeqLength];
    }
      
    // set of active users evolves according to the given probabilities
	for(uint iTime=_trainSeqLength+1;iTime<_frameLength;iTime++)
		for(uint iUser=0;iUser<_symbols.rows();iUser++)
		{   
			_usersActivity[iUser][iTime] = _usersActivityPdfs[iUser].sampleGivenItWas(_usersActivity[iUser][iTime-1]);             
			_symbols(iUser,_preambleLength+iTime) = _symbols(iUser,_preambleLength+iTime)*double(_usersActivity[iUser][iTime]);
			_isSymbolAccountedForDetection[iUser][iTime] = _usersActivity[iUser][iTime];
		}
            
#ifdef PRINT_INFO
    cout << "symbols after generating users activity" << endl << symbols << endl;
#endif
	
	do
	{	
		// spreading codes for the users are generated randomly
		_spreadingCodes = Util::sign(StatUtil::randnMatrix(_L,_N,0.0,1.0));

// 		// the spreading codes are normalized
// 		_spreadingCodes /= sqrt(_L);

// 		MatrixXd kasamiCodes (_L,_N);
// 		kasamiCodes <<	1,  -1,  -1, -1,   1,   1, 1,   1,   1, -1,  -1,   1, 1,   1,  -1, 1,  -1,   1, -1,  -1,  -1, -1,   1,  -1;
// 		_spreadingCodes = kasamiCodes;

#ifdef PRINT_CODES_INFO
		std::cout << "generated spreadingCodes..." << std::endl << _spreadingCodes << std::endl;
		std::cout << "are codes are orthogonal? " << Util::areColsOrthogonal(_spreadingCodes) << std::endl;
		std::cout << "are codes different and NOT opposite? " << Util::areColsDifferentAndNotOpposite(_spreadingCodes) << std::endl;
#endif
// 	} while(false);
	} while(!Util::areColsDifferentAndNotOpposite(_spreadingCodes));
	
	double thisChannelMinimumSIR;
	
	do
	{
		delete _channel;

// 		_channel = new MultiuserCDMAchannel(new ARchannel(_N,1,_m,_symbols.cols(),ARprocess(_powerProfile->generateChannelMatrix(_randomGenerator),ARcoefficients,ARvariance)),_spreadingCodes);
// 		channel = new MultiuserCDMAchannel(new TimeInvariantChannel(powerProfile->nInputs(),powerProfile->nOutputs(),m,symbols.cols(),MatrixXd::Ones(powerProfile->nOutputs(),powerProfile->nInputs())),_spreadingCodes);
		_channel = new MultiuserCDMAchannel(new BesselChannel(_N,1,_m,_symbols.cols(),_velocity,_carrierFrequency,_T,*_powerProfile),_spreadingCodes);
		
		// Signal to Interference Ratio's
		std::vector<double> SIRs = dynamic_cast<MultiuserCDMAchannel *> (_channel)->signalToInterferenceRatio(_iUserOfInterest);
		
		// the minimum is obtained
		thisChannelMinimumSIR = 10*log10(SIRs[std::min_element(SIRs.begin(),SIRs.end())-SIRs.begin()]);
		std::cout << "minimum SIR = " << thisChannelMinimumSIR << std::endl;

// 	} while(false);
	} while(thisChannelMinimumSIR<_minSignalToInterferenceRatio);	
	
	// the noise is generated
	_noise = new PowerProfileDependentNoise(_alphabet->variance(),_L,_channel->length(),*_powerProfile);
// 	_noise = new SingleUserChannelDependentNoise(_alphabet->variance(),_channel,_iUserOfInterest);
}

double CDMASystem::computeActivityDetectionErrorRate(MatrixXd sourceSymbols, MatrixXd detectedSymbols) const
{
	return computeSelectedUsersActivityDetectionErrorRate(sourceSymbols.row(_iUserOfInterest),detectedSymbols.row(_iUserOfInterest));
// 	return computeSelectedUsersActivityDetectionErrorRate(sourceSymbols,detectedSymbols);
}

double CDMASystem::computeSelectedUsersActivityDetectionErrorRate(MatrixXd sourceSymbols, MatrixXd detectedSymbols) const
{
	if (_symbolsDetectionWindowStart!=0)
		throw RuntimeException("CDMASystem::computeActivityDetectionErrorRate: this is only implemented when the starting point for SER computing is zero (the beginning of the frame).");

	// the dimensions of "sourceSymbols" are the same as those of "detectedSymbols"
	assert(sourceSymbols.rows() == detectedSymbols.rows() && sourceSymbols.cols()== detectedSymbols.cols());

	uint nSymbolsRows = detectedSymbols.rows();

	std::vector<std::vector<bool> > mask(nSymbolsRows,std::vector<bool>(_frameLength,true));

	// in order to compute the probability of activity detection it makes no difference the symbol detected: the only thing that matters is wether a symbol (any) was detected or not
	// for both the "sourceSymbols" and the "detectedSymbols", every symbol belonging to the alphabet is transformed into the "first" symbol of the alphabet
	for (uint i=0;i<sourceSymbols.rows();i++)
		for (uint j=0;j<sourceSymbols.cols();j++)
		{
			if (_alphabet->doesItBelong(sourceSymbols(i,j)))
				sourceSymbols(i,j) = _alphabet->operator[](0u);
			if (_alphabet->doesItBelong(detectedSymbols(i,j)))
				detectedSymbols(i,j) = _alphabet->operator[](0u);
		}

	double res = 0.0;

	for (uint iSignChange=1;iSignChange<_signChanges.size();iSignChange++)
	{
		res += (_signChanges[iSignChange]-_signChanges[iSignChange-1])*
				computeSERwithoutSolvingAmbiguity(
					sourceSymbols.block(0,_signChanges[iSignChange-1],nSymbolsRows,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
					Util::applyPermutationOnRows(detectedSymbols,_permutations[_piecesBestPermuationIndexes[iSignChange-1]],vector<int>(nSymbolsRows,+1)).block(0,_signChanges[iSignChange-1],nSymbolsRows,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
					Util::block(mask,0,_signChanges[iSignChange-1],nSymbolsRows,_signChanges[iSignChange]-_signChanges[iSignChange-1]));
	}

	res /= sourceSymbols.cols();

	return res;
}

void CDMASystem::beforeEndingAlgorithm()
{
	SMCSystem::beforeEndingAlgorithm();

	if(_algorithms[_iAlgorithm]->performsSymbolsDetection())
		_presentFramePeActivityDetection(_iSNR,_iAlgorithm) = computeActivityDetectionErrorRate(_symbols.block(0,_preambleLength,_N,_frameLength),_detectedSymbols);
	else
		_presentFramePeActivityDetection(_iSNR,_iAlgorithm) = -3.14;
	
	resetFramePieces();
}

void CDMASystem::onlyOnce()
{
	SMCSystem::onlyOnce();

	_presentFramePeActivityDetection = MatrixXd::Zero(_SNRs.size(),_algorithms.size());
	
	// the no sign changes is needed (the algorithm knows the channel), the value is set to 0
	_thisFrameNumberSignChanges = std::vector<std::vector<uint> >(_SNRs.size(),std::vector<uint>(_algorithms.size(),0u));
}

double CDMASystem::computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns)
{
	return computeSelectedUsersSER(sourceSymbols.row(_iUserOfInterest),detectedSymbols.row(_iUserOfInterest),Util::row(mask,_iUserOfInterest),iBestPermutation,bestPermutationSigns);	
// 	return computeSelectedUsersSER(sourceSymbols,detectedSymbols,mask,iBestPermutation,bestPermutationSigns);
}

double CDMASystem::computeSelectedUsersSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns)
{
	if (_symbolsDetectionWindowStart!=0)
		throw RuntimeException("CDMASystem::computeSER: this is only implemented when the starting point for SER computing is zero (the beginning of the frame).");

	uint nSymbolsRows = detectedSymbols.rows();
		
	// we dismiss the initialization in the constructor
	_piecesBestPermuationIndexes.clear();
	_piecesBestPermutationSigns.clear();

	double res = 0.0;

	// NOTE: even though ultimately only symbol vectors from "_symbolsDetectionWindowStart" will be taken into account for detection, we don't have to worry about it here, since the mask will take care of that. That being so, 
	// "_signChanges" actually contains all the sign changes that occur within the frame
	
	// if the algorithm performs channel estimation we use ITS channel estimations to split the frame and tackle the ambiguity problem
	if(_algorithms[_iAlgorithm]->performsChannelEstimation())
// 		_signChanges = Util::getZeroCrossings(_channel->getChannelMatrices(),_preambleLength,_frameLength);
// 		_signChanges = Util::getZeroCrossings(_algorithms[_iAlgorithm]->getEstimatedChannelMatrices(),_preambleLength,_frameLength);
		_signChanges = Util::getZeroCrossings(_algorithms[_iAlgorithm]->getEstimatedChannelMatrices(),_iUserOfInterest,_preambleLength,_frameLength);
	//...if it doesn't perform channel estimation, the ambiguity problem is gone
	else
	{
		_signChanges = std::vector<uint>(2);
		_signChanges[0] = 0; _signChanges[1] = _frameLength;
	}

	_thisFrameNumberSignChanges[_iSNR][_iAlgorithm] = _signChanges.size();
	
	uint overallNumberAccountedSymbols = 0;

	for (uint iSignChange=1;iSignChange<_signChanges.size();iSignChange++)
	{
		// we need to find out how many symbols are gonna be taken into account within this subframe
		uint thisSubframeNumberAccountedSymbols = 0;
		for (uint i=0;i<nSymbolsRows;i++)
			for (uint j=_signChanges[iSignChange-1];j<_signChanges[iSignChange];j++)
				thisSubframeNumberAccountedSymbols += mask[i][j];

		overallNumberAccountedSymbols += thisSubframeNumberAccountedSymbols;

		res += thisSubframeNumberAccountedSymbols*
				BaseSystem::computeSER(sourceSymbols.block(0,_signChanges[iSignChange-1],nSymbolsRows,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
										detectedSymbols.block(0,_signChanges[iSignChange-1],nSymbolsRows,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
										Util::block(mask,0,_signChanges[iSignChange-1],nSymbolsRows,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
										iBestPermutation,bestPermutationSigns);
		
#ifdef DEBUG
		cout << "iBestPermutation = " << iBestPermutation << endl << "bestPermutationSigns:" << endl << bestPermutationSigns << endl;
		cout << detectedSymbols.block(0,_frameLength-50,nSymbolsRows,50);
#endif
		assert(iBestPermutation==0);
		
		// we need to store which the best permutations were along with their corresponding signs since it will be needed later by "computeActivityDetectionErrorRate" and "computeMSE"
		_piecesBestPermuationIndexes.push_back(iBestPermutation);
		_piecesBestPermutationSigns.push_back(bestPermutationSigns);
  }

	res /= overallNumberAccountedSymbols;

	return res;
}

double CDMASystem::computeMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices) const
{
	return computeSelectedUsersMSE(Util::keepCol(realChannelMatrices,_iUserOfInterest),Util::keepCol(estimatedChannelMatrices,_iUserOfInterest));
// 	return computeSelectedUsersMSE(realChannelMatrices,estimatedChannelMatrices);
}


double CDMASystem::computeSelectedUsersMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices) const
{
  if(realChannelMatrices.size()!=estimatedChannelMatrices.size())
	throw RuntimeException("CDMASystem::computeMSE: different number of real channel matrices than estimated channel matrices.");

  // we know that the received channel matrices go from time instant preambleLength+MSEwindowStart (both known parameters of the system) to the end of the frame. The first matrix of realChannelMatrices/estimatedChannelMatrices corresponds to the time instant:
  uint iChannelMatricesStart = _preambleLength+_MSEwindowStart;
  
  // we find which interval "iChannelMatricesStart" belongs to ("_signChanges" was computed previously in "computeSER")
  uint iSignChange = 1;
  while(_signChanges[iSignChange]<=iChannelMatricesStart && iSignChange<_signChanges.size())
	iSignChange++;
  

  std::vector<MatrixXd> toCheckRealChannelMatrices(realChannelMatrices.begin(),realChannelMatrices.begin()+_signChanges[iSignChange]-iChannelMatricesStart);
  std::vector<MatrixXd> toCheckEstimatedChannelMatrices(estimatedChannelMatrices.begin(),estimatedChannelMatrices.begin()+_signChanges[iSignChange]-iChannelMatricesStart);
  
  double res = (_signChanges[iSignChange]-iChannelMatricesStart)*BaseSystem::computeMSE(toCheckRealChannelMatrices,toCheckEstimatedChannelMatrices,
								_permutations[_piecesBestPermuationIndexes[iSignChange-1]],_piecesBestPermutationSigns[iSignChange-1]);
  iSignChange++;
  
  for(;iSignChange<_signChanges.size();iSignChange++)
  {
	toCheckRealChannelMatrices = std::vector<MatrixXd>(realChannelMatrices.begin()+_signChanges[iSignChange-1]-iChannelMatricesStart,realChannelMatrices.begin()+_signChanges[iSignChange]-iChannelMatricesStart);
	toCheckEstimatedChannelMatrices = std::vector<MatrixXd>(estimatedChannelMatrices.begin()+_signChanges[iSignChange-1]-iChannelMatricesStart,estimatedChannelMatrices.begin()+_signChanges[iSignChange]-iChannelMatricesStart);
	
	res += (_signChanges[iSignChange]-iChannelMatricesStart)*BaseSystem::computeMSE(toCheckRealChannelMatrices,toCheckEstimatedChannelMatrices,
								_permutations[_piecesBestPermuationIndexes[iSignChange-1]],_piecesBestPermutationSigns[iSignChange-1]);
  }

  res /= realChannelMatrices.size();

  return res;
}

void CDMASystem::resetFramePieces()
{
	// by default (when "computeSER" is not called), we assume that there is no sign changes (the frame is not split)...
	_signChanges = std::vector<uint>(2);
	_signChanges[0] = 0; _signChanges[1] = _frameLength;
	
	// ...the best permutation is the first one...
	_piecesBestPermuationIndexes = std::vector<uint>(1,0);
	
	// ...and its corresponding signs are all +1
	_piecesBestPermutationSigns = std::vector<std::vector<int> >(1,std::vector<int>(_permutations[0].size(),+1));
}

void CDMASystem::storeFrameResults()
{
    SMCSystem::storeFrameResults();

    _peActivityDetectionFrames.push_back(_presentFramePeActivityDetection);
	_everyFrameUsersActivity.push_back(_usersActivity);
	_everyFrameSpreadingCodes.push_back(_spreadingCodes);
	_everyFrameNumberSignChanges.push_back(_thisFrameNumberSignChanges);
}

void CDMASystem::saveFrameResults()
{
    SMCSystem::saveFrameResults();
    Octave::eigenToOctaveFileStream(_peActivityDetectionFrames,"peActivityDetectionFrames",_f);
    Octave::eigenToOctaveFileStream(_spreadingCodes,"spreadingCodes",_f);
	Octave::toOctaveFileStream(_nSurvivors,"nSurvivors",_f);
	Octave::toOctaveFileStream(_userPersistenceProb,"userPersistenceProb",_f);
	Octave::toOctaveFileStream(_newActiveUserProb,"newActiveUserProb",_f);
	Octave::toOctaveFileStream(_userPriorProb,"userPriorProb",_f);
	Octave::toOctaveFileStream(_everyFrameUsersActivity,"usersActivity",_f);
	Octave::toOctaveFileStream(_signChanges,"signChanges",_f);
	Octave::toOctaveFileStream(_minSignalToInterferenceRatio,"minSignalToInterferenceRatio",_f);
	Octave::eigenToOctaveFileStream(_everyFrameSpreadingCodes,"everyFrameSpreadingCodes",_f);
	Octave::toOctaveFileStream(_everyFrameNumberSignChanges,"everyFrameNumberSignChanges",_f);
	
	Octave::toOctaveFileStream(_velocity,"velocity",_f);
	Octave::toOctaveFileStream(_carrierFrequency,"carrierFrequency",_f);
	Octave::toOctaveFileStream(_symbolRate,"symbolRate",_f);
	Octave::toOctaveFileStream(_T,"T",_f);
	
// 	ARprocess miar(StatUtil::randnMatrix(1,3,0.0,1.0),2,_velocity,_carrierFrequency,_T);
// 	std::vector<MatrixXd> m;
// 	
// 	for(uint i=0;i<1000;i++)
// 		m.push_back(miar.nextMatrix());
// 	
// 	Octave::eigenToOctaveFileStream(m,"matrices",_f);
}
