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
#include <CDMAunknownActiveUsersSISopt.h>
#include <TimeInvariantChannel.h>
#include <MultiuserCDMAchannel.h>
#include <ViterbiAlgorithmWithAprioriProbabilities.h>
#include <PSPAlgorithmWithAprioriProbabilities.h>

#include <math.h>

#include <bashcolors.h>

#define PRINT_CODES_INFO
// #define PRINT_ACTIVITY_SAMPLING

// #define PRINT_INFO

// #define DEBUG_SER
// #define DEBUG_MSE
// #define DEBUG_ACTIVITY_DETECTION_ERROR_RATE

CDMASystem::CDMASystem(): SMCSystem()
,_piecesInfoAvailable(false)
,_userPersistenceProb(0.99),_newActiveUserProb(0.01),_userPriorProb(0.5)
// ,userPersistenceProb(0.8),newActiveUserProb(0.2),userPriorProb(1.0)
// ,userPersistenceProb(1.0),newActiveUserProb(0.2),userPriorProb(1.0)
,_usersActivityPdfs(_N,UsersActivityDistribution(_userPersistenceProb,_newActiveUserProb,_userPriorProb))
// ,maximumRatioThresholdInDBs(15)
,_maximumRatioThresholdInDBs(20)
{
    if (_m!=1)
        throw RuntimeException("CDMASystem::CDMASystem: channel is not flat.");

	// first user starts transmitting something
	_usersActivityPdfs[0].setApriori(1.0);
	
    // spreading spreadingCodes for the users are generated randomly
//     _spreadingCodes = StatUtil::randnMatrix(L,N,0.0,1.0);
//     _spreadingCodes = Util::sign(_spreadingCodes);

	MatrixXd kasamiCodes (_L,_N);
	
	kasamiCodes <<	 1,  -1,  -1,
					-1,   1,   1,
					 1,   1,   1,
					-1,  -1,   1,
					 1,   1,  -1,
					 1,  -1,   1,
					-1,  -1,  -1,
					-1,   1,  -1;

	_spreadingCodes = kasamiCodes;
	
	// the spreading codes are normalized
// 	_spreadingCodes /= sqrt(L);
    
#ifdef PRINT_CODES_INFO
    cout << "generated spreadingCodes..." << endl << _spreadingCodes << endl;
	cout << "are codes are ok? " << areSequencesOrthogonal(_spreadingCodes) << endl;
#endif

  _nSurvivors = 2;
//   _nSurvivors = 8;
//   _nSurvivors = 10;
//   _nSurvivors = 20;
//   _nSurvivors = 40;

    // AR process parameters
    ARcoefficients = vector<double>(2);
    ARcoefficients[0] = 0.59999;
    ARcoefficients[1] = 0.39999;
    ARvariance=0.0001;    
    
//     // AR process parameters
//     ARcoefficients = vector<double>(1);
//     ARcoefficients[0] = 0.99999;
//     ARvariance=0.0001;
    
    // a flat power profile is generated. Notice:
    //      i) that m should be 1, otherwise an exception would have been thrown
    //     ii) we only need to generate a coefficient per user, i.e., a 1xN vector
    _powerProfile = new FlatPowerProfile(1,_N,_m,1.0);
    
    _cdmaKalmanEstimator = new CDMAKalmanEstimator(_powerProfile->means(),_powerProfile->variances(),ARcoefficients,ARvariance,_spreadingCodes);
    _cdmaKnownChannelChannelMatrixEstimator = NULL;
	
    _mmseDetector = new MMSEDetector(_L,_N,_alphabet->variance(),_N);

	// bessel channel parameters
    _velocity = 180/3.6; // (m/s)
    _carrierFrequency = 2e9; // (Hz)
    _symbolRate = 500e3; // (Hz)

    _T = 1.0/_symbolRate; // (s)
	
	_maxCoefficientsRatiosInDBs.reserve(_nFrames);
	
	_peActivityDetectionFrames.reserve(_nFrames);
	
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
	  nParticles = (int)pow((double)_alphabet->length()+1,_N)*_nSurvivors;
	  
	  cout << " to " << nParticles << endl;
    }

    if(_adjustSurvivorsFromParticlesNumber)
    {
	  cout << "Number of survivors adjusted from " << _nSurvivors;
	  _nSurvivors = int(ceil(double(nParticles)/pow((double)_alphabet->length()+1,double(_N))));
	  cout << " to " << _nSurvivors << endl;
    }
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
    _algorithms.push_back(new KnownFlatChannelOptimalAlgorithm ("CDMA optimal with known channel BUT no knowledge of users activity probabilities",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,*_channel,_preambleLength));
    
    _algorithms.push_back(new KnownFlatChannelAndActiveUsersOptimalAlgorithm ("CDMA optimal (known channel and active users)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,*_channel,_preambleLength,_usersActivity));    
       
    // the channel is different in each frame, so the estimator that knows the channel must be rebuilt every frame
    delete _cdmaKnownChannelChannelMatrixEstimator;
    _cdmaKnownChannelChannelMatrixEstimator = new CDMAKnownChannelChannelMatrixEstimator(_channel,_preambleLength,_N,_spreadingCodes);
     
//     algorithms.push_back(new CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge ("CDMA SIS-opt with no knowledge of users activity pdf",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances()));

    _algorithms.push_back(new CDMAunknownActiveUsersSISopt ("CDMA SIS-opt",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_preamble,_d,nParticles,algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_usersActivityPdfs));

    _algorithms.push_back(new UnknownActiveUsersLinearFilterBasedSMCAlgorithm ("CDMA SIS Linear Filters",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_mmseDetector,_preamble,_d,nParticles,algoritmoRemuestreo,_powerProfile->means(),_powerProfile->variances(),_usersActivityPdfs));

	_algorithms.push_back(new ViterbiAlgorithmWithAprioriProbabilities("Viterbi with a priori probabilities (known channel)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (_channel)),_preamble,_d,_usersActivityPdfs));
	
	_algorithms.push_back(new PSPAlgorithmWithAprioriProbabilities("PSP",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_preamble,_d,_iLastSymbolVectorToBeDetected+_d,_nSurvivors,_usersActivityPdfs));
	
	_algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*_alphabet,_L,1,_N,_iLastSymbolVectorToBeDetected,_m,_cdmaKalmanEstimator,_preamble,_symbols));
}

void CDMASystem::beforeEndingFrame()
{
    SMCSystem::beforeEndingFrame();

	// the maximum ratio of this frame is added to the vector
	_maxCoefficientsRatiosInDBs.push_back(_maximumRatio);
	
    _peActivityDetectionFrames.push_back(_presentFramePeActivityDetection);
    Util::matricesVectorToOctaveFileStream(_peActivityDetectionFrames,"peActivityDetectionFrames",_f);
	
	Util::scalarsVectorToOctaveFileStream(_maxCoefficientsRatiosInDBs,"maxCoefficientsRatiosInDBs",_f);
    Util::matrixToOctaveFileStream(_spreadingCodes,"spreadingCodes",_f);
	Util::scalarToOctaveFileStream(_maximumRatioThresholdInDBs,"maximumRatioThresholdInDBs",_f);
	Util::scalarToOctaveFileStream(_nSurvivors,"nSurvivors",_f);
	Util::scalarToOctaveFileStream(_userPersistenceProb,"userPersistenceProb",_f);
	Util::scalarToOctaveFileStream(_newActiveUserProb,"newActiveUserProb",_f);
	Util::scalarToOctaveFileStream(_userPriorProb,"userPriorProb",_f);
}

void CDMASystem::buildChannel()
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
	for(int iTime=_trainSeqLength+1;iTime<_frameLength;iTime++)
		for(int iUser=0;iUser<_symbols.rows();iUser++)
		{   
			_usersActivity[iUser][iTime] = _usersActivityPdfs[iUser].sampleGivenItWas(_usersActivity[iUser][iTime-1]);             
			_symbols(iUser,_preambleLength+iTime) = _symbols(iUser,_preambleLength+iTime)*double(_usersActivity[iUser][iTime]);
			_isSymbolAccountedForDetection[iUser][iTime] = _usersActivity[iUser][iTime];
		}
            
#ifdef PRINT_INFO
    cout << "symbols after generating users activity" << endl << symbols << endl;
#endif    

	bool coefficientsSignChangeHappened;
	std::vector<uint> channelCoefficientsSignChanges;
	
	do
	{
		delete _channel;

// 		_channel = new MultiuserCDMAchannel(new ARchannel(_N,1,_m,_symbols.cols(),ARprocess(_powerProfile->generateChannelMatrix(_randomGenerator),ARcoefficients,ARvariance)),_spreadingCodes);
// 		channel = new MultiuserCDMAchannel(new TimeInvariantChannel(powerProfile->nInputs(),powerProfile->nOutputs(),m,symbols.cols(),MatrixXd::Ones(powerProfile->nOutputs(),powerProfile->nInputs())),_spreadingCodes);
		_channel = new MultiuserCDMAchannel(new BesselChannel(_N,1,_m,_symbols.cols(),_velocity,_carrierFrequency,_T,*_powerProfile),_spreadingCodes);
	  
		channelCoefficientsSignChanges = _channel->getInputsZeroCrossings(_preambleLength,_frameLength);
		coefficientsSignChangeHappened = channelCoefficientsSignChanges.size() > 2;
		if(coefficientsSignChangeHappened)
		{
			cout << "any of the channel coefficients changed sign at " << channelCoefficientsSignChanges[1] << endl;
// 			cout << "channel" << endl << _channel->at(_preambleLength) << endl;
// 			getchar();
		}

	} while(!isChannelOk(_channel));
// 	} while(!isChannelOk(_channel) || !coefficientsSignChangeHappened); // los coeficientes del canal cambian de signo
// 	} while(!isChannelOk(_channel) || coefficientsSignChangeHappened); // los coeficientes del canal NO cambian de signo
}

bool CDMASystem::areSequencesOrthogonal(const MatrixXd &spreadingCodes)
{
	uint L = spreadingCodes.rows();
	uint nCodes = spreadingCodes.cols();

	for (uint iOneCode=0;iOneCode<nCodes;iOneCode++)
		for (uint iOtherCode=iOneCode+1;iOtherCode<nCodes;iOtherCode++)
		{
			int sum = 0;
			
			for (uint i=0;i<L;i++)
				for (uint j=0;j<L;j++)
					sum += spreadingCodes(i,iOneCode)*spreadingCodes(j,iOtherCode);
			
			if (sum!=0)
				return false;
		}

	return true;
}

double CDMASystem::computeActivityDetectionErrorRate(MatrixXd sourceSymbols, MatrixXd detectedSymbols) const
{
	if(!_piecesInfoAvailable)
		throw RuntimeException("CDMASystem::computeMSE: pieces information is not available.");

	if (_symbolsDetectionWindowStart!=0)
		throw RuntimeException("CDMASystem::computeActivityDetectionErrorRate: this is only implemented when the starting point for SER computing is zero (the beginning of the frame).");

	if (detectedSymbols.rows()==0)
		return -1.0;

	if (sourceSymbols.rows()!= detectedSymbols.rows())
	{
		cout << "sourceSymbols.rows() = " << sourceSymbols.rows() << " detectedSymbols.rows() = " << detectedSymbols.rows() << endl;
		throw RuntimeException("CDMASystem::computeActivityDetectionER: matrix row numbers differ.");
	}

	if (sourceSymbols.cols()!= detectedSymbols.cols())
	{
		cout << "sourceSymbols.cols() = " << sourceSymbols.cols() << " detectedSymbols.cols() = " << detectedSymbols.cols() << endl;
		throw RuntimeException("CDMASystem::computeActivityDetectionER: matrix column numbers differ.");
	}

	vector<vector<bool> > mask(_N,vector<bool>(_frameLength,true));

	for (int iTime=0;iTime<_symbolsDetectionWindowStart;iTime++)
		for (uint iInput=0;iInput<_N;iInput++)
			mask[iInput][iTime] = false;

	// in order to compute the probability of activity detection it makes no difference the symbol detected: the only thing that matters is wether a symbol (any) was detected or not
	// for both the "sourceSymbols" and the "detectedSymbols", every symbol belonging to the alphabet is transformed into the "first" symbol of the alphabet
	for (int i=0;i<sourceSymbols.rows();i++)
		for (int j=0;j<sourceSymbols.cols();j++)
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
					sourceSymbols.block(0,_signChanges[iSignChange-1],_N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
					Util::applyPermutationOnRows(detectedSymbols,_permutations[_piecesBestPermuationIndexes[iSignChange-1]],vector<int>(_N,1)).block(0,_signChanges[iSignChange-1],_N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
					Util::block(mask,0,_signChanges[iSignChange-1],_N,_signChanges[iSignChange]-_signChanges[iSignChange-1]));
	}

	res /= sourceSymbols.cols();

	return res;
}

void CDMASystem::beforeEndingAlgorithm()
{
	SMCSystem::beforeEndingAlgorithm();

	// FIXME: "computeActivityDetectionErrorRate" depends on a result obtained during the execution of "computeSER", which is run by BaseSystem :-S
	double peActivity = computeActivityDetectionErrorRate(_symbols.block(0,_preambleLength,_N,_frameLength),_detectedSymbols);

	// the activity detection error probability is accumulated
	_presentFramePeActivityDetection(_iSNR,_iAlgorithm) = peActivity;
}

void CDMASystem::onlyOnce()
{
	SMCSystem::onlyOnce();

	_presentFramePeActivityDetection = MatrixXd::Zero(_SNRs.size(),_algorithms.size());
}

bool CDMASystem::isChannelOk(const MIMOChannel * const channel)
{
	double thisChannelMatrixMaximumRatio;
	int iChannelMatrix;

	// we check if the channel is really bad (severe near-far issues)...
	_maximumRatio = 20*log10(Util::maxCoefficientsRatio(channel->at(_preambleLength)));

	// all the channel matrices contained in this channel are checked
	for (iChannelMatrix=_preambleLength+1;iChannelMatrix<channel->length();iChannelMatrix++)
	{
		// check if how large is the difference of power between coefficients
		thisChannelMatrixMaximumRatio = 20*log10(Util::maxCoefficientsRatio(channel->at(iChannelMatrix)));
		if (thisChannelMatrixMaximumRatio < _maximumRatio)
			_maximumRatio = thisChannelMatrixMaximumRatio;
	}

	cout << COLOR_WHITE << "the max difference among coefficients in dBs: " << COLOR_NORMAL << _maximumRatio << endl;

	if (_maximumRatio > _maximumRatioThresholdInDBs)
	{
		cout << COLOR_PINK << "the max difference among coefficients in dBs is " << COLOR_NORMAL << _maximumRatio << COLOR_PINK << "...channel is NOT ok!!" << COLOR_NORMAL << endl;
		return false;
	}

	return true;
}

double CDMASystem::computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns)
{
	if (_symbolsDetectionWindowStart!=0)
		throw RuntimeException("CDMASystem::computeSER: this is only implemented when the starting point for SER computing is zero (the beginning of the frame).");

	if (detectedSymbols.rows()==0)
	{
		// this is used later by "computeActivityDetectionErrorRate" and "computeMSE"
		_piecesBestPermuationIndexes = std::vector<uint>(_signChanges.size()-1,0);
		_piecesBestPermutationSigns = std::vector<std::vector<int> >(_signChanges.size()-1,std::vector<int>(_N,1));
		return -1.0;
	}

	_piecesBestPermuationIndexes.clear();
	_piecesBestPermutationSigns.clear();

	double res = 0.0;

// 	_signChanges = channel->getInputsZeroCrossings(preambleLength+symbolsDetectionWindowStart,frameLength-symbolsDetectionWindowStart);
	//								^
	//								|
	// even though ultimately only symbol vectors from "symbolsDetectionWindowStart" will be taken into account for detection, we don't
	// have to worry about it here, since the mask will take care of that. That being so, "_signChanges" actually contains all the sign changes
	// that occur within the frame
	_signChanges = _channel->getInputsZeroCrossings(_preambleLength,_frameLength);

	uint overallNumberAccountedSymbols = 0;

	for (uint iSignChange=1;iSignChange<_signChanges.size();iSignChange++)
	{
		// we need to find out how many symbols are gonna be taken into account within this subframe
		uint thisSubframeNumberAccountedSymbols = 0;
		for (uint i=0;i<_N;i++)
			for (uint j=_signChanges[iSignChange-1];j<_signChanges[iSignChange];j++)
				thisSubframeNumberAccountedSymbols += mask[i][j];

		overallNumberAccountedSymbols += thisSubframeNumberAccountedSymbols;

		res += thisSubframeNumberAccountedSymbols*
				BaseSystem::computeSER(sourceSymbols.block(0,_signChanges[iSignChange-1],_N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
										detectedSymbols.block(0,_signChanges[iSignChange-1],_N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
										Util::block(mask,0,_signChanges[iSignChange-1],_N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
										iBestPermutation,bestPermutationSigns);
		
		// we need to store which the best permutations were along with their corresponding signs since it will be needed later by "computeActivityDetectionErrorRate" and "computeMSE"
		_piecesBestPermuationIndexes.push_back(iBestPermutation);
		_piecesBestPermutationSigns.push_back(bestPermutationSigns);
  }

	res /= overallNumberAccountedSymbols;

	// the information concerning how the data frame must be cut to measure performance becomes available
	_piecesInfoAvailable = true;
	
	return res;
}

double CDMASystem::computeMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices) const
{
	if(!_piecesInfoAvailable)
		throw RuntimeException("CDMASystem::computeMSE: pieces information is not available.");

  if(estimatedChannelMatrices.size()==0)
	return -1.0;
  
  if(realChannelMatrices.size()!=estimatedChannelMatrices.size())
	throw RuntimeException("CDMASystem::computeMSE: different number of real channel matrices than estimated channel matrices.");

  // we know that the received channel matrices go from time instant preambleLength+MSEwindowStart (both known parameters of the system) 
  // to the end of the frame. The first matrix of realChannelMatrices/estimatedChannelMatrices corresponds to the time instant:
  uint iChannelMatricesStart = _preambleLength+_MSEwindowStart;
  
  // we find which interval "iChannelMatricesStart" belongs to ("_signChanges" was computed previously in "computeSER")
  uint iSignChange = 1;
  while(_signChanges[iSignChange]<=iChannelMatricesStart && iSignChange<_signChanges.size())
	iSignChange++;
  
#ifdef DEBUG_MSE
  cout << "computing MSE between " << iChannelMatricesStart << " and " << _signChanges[iSignChange] << " ( " << (_signChanges[iSignChange]-iChannelMatricesStart) << " matrices)" << endl;
  cout << "the permutation is " << _piecesBestPermuationIndexes[iSignChange-1] << endl;
  cout << "and the signs" << endl;
  Util::print(_piecesBestPermutationSigns[iSignChange-1]);
  cout << endl;
#endif

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

#ifdef DEBUG_MSE
	cout << "computing MSE between " << _signChanges[iSignChange-1] << " and " << _signChanges[iSignChange] << " ( " << (_signChanges[iSignChange]-_signChanges[iSignChange-1]) << " matrices)" << endl;
#endif
  }

  res /= realChannelMatrices.size();

  return res;
}