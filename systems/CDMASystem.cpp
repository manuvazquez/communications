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
,_userPersistenceProb(0.99),_newActiveUserProb(0.01),_userPriorProb(0.5)
// ,userPersistenceProb(0.8),newActiveUserProb(0.2),userPriorProb(1.0)
// ,userPersistenceProb(1.0),newActiveUserProb(0.2),userPriorProb(1.0)
,_usersActivityPdfs(N,UsersActivityDistribution(_userPersistenceProb,_newActiveUserProb,_userPriorProb))
// ,maximumRatioThresholdInDBs(15)
,_maximumRatioThresholdInDBs(20)
// ,_extendedAlphabet(alphabet->buildNewAlphabetByAddingSymbol(0.0))
{
    if(m!=1)
        throw RuntimeException("CDMASystem::CDMASystem: channel is not flat.");

	// first users starts transmitting something
	_usersActivityPdfs[0].setApriori(1.0);
	
    // spreading spreadingCodes for the users are generated randomly
//     _spreadingCodes = StatUtil::randnMatrix(L,N,0.0,1.0);
//     _spreadingCodes = Util::sign(_spreadingCodes);

	MatrixXd kasamiCodes (L,N);
	
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

//   nSurvivors = 2;
//   nSurvivors = 8;
//   nSurvivors = 10;
//   nSurvivors = 20;
  _nSurvivors = 40;

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
    powerProfile = new FlatPowerProfile(1,N,m,1.0);
    
    _cdmaKalmanEstimator = new CDMAKalmanEstimator(powerProfile->means(),powerProfile->variances(),ARcoefficients,ARvariance,_spreadingCodes);
    _cdmaKnownChannelChannelMatrixEstimator = NULL;
	
    _mmseDetector = new MMSEDetector(L,N,alphabet->variance(),N);

	// bessel channel parameters
    _velocity = 50/3.6; // (m/s)
    _carrierFrequency = 2e9; // (Hz)
    _symbolRate = 500e3; // (Hz)

    _T = 1.0/_symbolRate; // (s)
	
	_maxCoefficientsRatiosInDBs.reserve(nFrames);
	
	_peActivityDetectionFrames.reserve(nFrames);
	
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
	  nParticles = (int)pow((double)alphabet->length()+1,N)*_nSurvivors;
	  
	  cout << " to " << nParticles << endl;
    }

    if(_adjustSurvivorsFromParticlesNumber)
    {
	  cout << "Number of survivors adjusted from " << _nSurvivors;
	  _nSurvivors = int(ceil(double(nParticles)/pow((double)alphabet->length()+1,double(N))));
	  cout << " to " << _nSurvivors << endl;
    }
}


CDMASystem::~CDMASystem()
{
    delete powerProfile;
    delete _cdmaKalmanEstimator;
    delete _cdmaKnownChannelChannelMatrixEstimator;
    delete _mmseDetector;
}

void CDMASystem::AddAlgorithms()
{
    algorithms.push_back(new KnownFlatChannelOptimalAlgorithm ("CDMA optimal with known channel BUT no knowledge of users activity probabilities",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,*channel,preambleLength));
    
    algorithms.push_back(new KnownFlatChannelAndActiveUsersOptimalAlgorithm ("CDMA optimal (known channel and active users)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,*channel,preambleLength,_usersActivity));    
       
    // the channel is different in each frame, so the estimator that knows the channel must be rebuilt every frame
    delete _cdmaKnownChannelChannelMatrixEstimator;
    _cdmaKnownChannelChannelMatrixEstimator = new CDMAKnownChannelChannelMatrixEstimator(channel,preambleLength,N,_spreadingCodes);
     
// //     algorithms.push_back(new CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge ("CDMA SIS-opt with no knowledge of users activity pdf",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances()));

    algorithms.push_back(new CDMAunknownActiveUsersSISopt ("CDMA SIS-opt",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,_cdmaKalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),_usersActivityPdfs));

    algorithms.push_back(new UnknownActiveUsersLinearFilterBasedSMCAlgorithm ("CDMA SIS Linear Filters",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,_cdmaKalmanEstimator,_mmseDetector,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),_usersActivityPdfs));

	algorithms.push_back(new ViterbiAlgorithmWithAprioriProbabilities("Viterbi (known channel)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d,_usersActivityPdfs));
	
	algorithms.push_back(new PSPAlgorithmWithAprioriProbabilities("PSP",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,_cdmaKalmanEstimator,preamble,d,iLastSymbolVectorToBeDetected+d,_nSurvivors,_usersActivityPdfs));
	
	algorithms.push_back(new KnownSymbolsKalmanBasedChannelEstimatorAlgorithm("Kalman Filter (Known Symbols)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,_cdmaKalmanEstimator,preamble,symbols));
}

void CDMASystem::BeforeEndingFrame()
{
    SMCSystem::BeforeEndingFrame();

	// the maximum ratio of this frame is added to the vector
	_maxCoefficientsRatiosInDBs.push_back(_maximumRatio);
	
    _peActivityDetectionFrames.push_back(_presentFramePeActivityDetection);
    Util::matricesVectorToOctaveFileStream(_peActivityDetectionFrames,"peActivityDetectionFrames",f);
	
	Util::scalarsVectorToOctaveFileStream(_maxCoefficientsRatiosInDBs,"maxCoefficientsRatiosInDBs",f);
    Util::matrixToOctaveFileStream(_spreadingCodes,"spreadingCodes",f);
	Util::scalarToOctaveFileStream(_maximumRatioThresholdInDBs,"maximumRatioThresholdInDBs",f);
	Util::scalarToOctaveFileStream(_nSurvivors,"nSurvivors",f);
	Util::scalarToOctaveFileStream(_userPersistenceProb,"userPersistenceProb",f);
	Util::scalarToOctaveFileStream(_newActiveUserProb,"newActiveUserProb",f);
	Util::scalarToOctaveFileStream(_userPriorProb,"userPriorProb",f);
}

void CDMASystem::BuildChannel()
{
    // when users are not transmitting, their symbols are zero
    _usersActivity = vector<vector<bool> >(symbols.rows(),vector<bool>(frameLength));
    
    // at the first time instant the prior probability is used to decide which users are active
    for(uint iUser=0;iUser<static_cast<uint>(symbols.rows());iUser++)
    {
        _usersActivity[iUser][trainSeqLength] = _usersActivityPdfs[iUser].sampleFromPrior();
#ifdef PRINT_ACTIVITY_SAMPLING
		cout << "user " << iUser << ": " << _usersActivity[iUser][trainSeqLength] << endl;
#endif
        symbols(iUser,preambleLength+trainSeqLength) = double(_usersActivity[iUser][trainSeqLength])*symbols(iUser,preambleLength+trainSeqLength);
        isSymbolAccountedForDetection[iUser][trainSeqLength] = _usersActivity[iUser][trainSeqLength];
    }
      
    // set of active users evolves according to the given probabilities
    for(int iTime=trainSeqLength+1;iTime<frameLength;iTime++)    
        for(int iUser=0;iUser<symbols.rows();iUser++)
        {   
            _usersActivity[iUser][iTime] = _usersActivityPdfs[iUser].sampleGivenItWas(_usersActivity[iUser][iTime-1]);             
            symbols(iUser,preambleLength+iTime) = symbols(iUser,preambleLength+iTime)*double(_usersActivity[iUser][iTime]);
            isSymbolAccountedForDetection[iUser][iTime] = _usersActivity[iUser][iTime];
        }
            
#ifdef PRINT_INFO
    cout << "symbols after generating users activity" << endl << symbols << endl;
#endif    

	bool aCoefficientChangesSign;

	do
	{
	  delete channel;

// 	  channel = new MultiuserCDMAchannel(new ARchannel(N,1,m,symbols.cols(),ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance)),_spreadingCodes);
	  
// 	  channel = new MultiuserCDMAchannel(new TimeInvariantChannel(powerProfile->nInputs(),powerProfile->nOutputs(),m,symbols.cols(),MatrixXd::Ones(powerProfile->nOutputs(),powerProfile->nInputs())),_spreadingCodes);
	  
	  channel = new MultiuserCDMAchannel(new BesselChannel(N,1,m,symbols.cols(),_velocity,_carrierFrequency,_T,*powerProfile),_spreadingCodes);
	  
	  aCoefficientChangesSign = channel->getInputsZeroCrossings(preambleLength,frameLength).size() > 2;
// 	  cout << "aCoefficientChangesSign = " << aCoefficientChangesSign << endl;

// 	} while(!isChannelOk(channel) || !aCoefficientChangesSign); // los coeficientes del canal cambian de signo
	} while(!isChannelOk(channel) || aCoefficientChangesSign); // los coeficientes del canal NO cambian de signo
}

bool CDMASystem::areSequencesOrthogonal(const MatrixXd &spreadingCodes)
{
  int L = spreadingCodes.rows();
  int nCodes = spreadingCodes.cols();
  
  for(int iOneCode=0;iOneCode<nCodes;iOneCode++)
	for(int iOtherCode=iOneCode+1;iOtherCode<nCodes;iOtherCode++)
	{
	  int sum = 0;
	  
	  for(int i=0;i<L;i++)
		for(int j=0;j<L;j++)
		  sum += spreadingCodes(i,iOneCode)*spreadingCodes(j,iOtherCode);
		
	  if(sum!=0)
		return false;
	}
	
  return true;
}

double CDMASystem::computeActivityDetectionErrorRate(MatrixXd sourceSymbols, MatrixXd detectedSymbols) const
{
  if(symbolsDetectionWindowStart!=0)
	throw RuntimeException("CDMASystem::computeActivityDetectionErrorRate: this is only implemented when the starting point for SER computing is zero (the beginning of the frame).");

  if(detectedSymbols.rows()==0)
	return -1.0;
  
#ifdef DEBUG
	cout << "source symbols" << endl << sourceSymbols << endl;
	cout << "detected symbols" << endl << detectedSymbols << endl;
#endif

  if(sourceSymbols.rows()!= detectedSymbols.rows())
  {
	  cout << "sourceSymbols.rows() = " << sourceSymbols.rows() << " detectedSymbols.rows() = " << detectedSymbols.rows() << endl;
	  throw RuntimeException("CDMASystem::computeActivityDetectionER: matrix row numbers differ.");
  }

  if(sourceSymbols.cols()!= detectedSymbols.cols())
  {
	cout << "sourceSymbols.cols() = " << sourceSymbols.cols() << " detectedSymbols.cols() = " << detectedSymbols.cols() << endl;    
	throw RuntimeException("CDMASystem::computeActivityDetectionER: matrix column numbers differ.");
  }

  vector<vector<bool> > mask(N,vector<bool>(frameLength,true));
  
//   return computeSERwithoutSolvingAmbiguity(sourceSymbols,detectedSymbols,mask);

  for(int iTime=0;iTime<symbolsDetectionWindowStart;iTime++)
	  for(int iInput=0;iInput<N;iInput++)
		  mask[iInput][iTime] = false;        

  // in order to compute the probability of activity detection it makes no difference the symbol detected
  // the only thing that matters is wether a symbol (any) was detected or not
  for(int i=0;i<sourceSymbols.rows();i++)
	for(int j=0;j<sourceSymbols.cols();j++)
	{
	  if(alphabet->doesItBelong(sourceSymbols(i,j)))
		sourceSymbols(i,j) = alphabet->operator[](0);
	  if(alphabet->doesItBelong(detectedSymbols(i,j)))
		detectedSymbols(i,j) = alphabet->operator[](0);
	}

#ifdef DEBUG
	cout << "despues" << endl;
	cout << "source symbols" << endl << sourceSymbols << endl;
	cout << "detected symbols" << endl << detectedSymbols << endl;
#endif

  double res = 0.0;

  for(uint iSignChange=1;iSignChange<_signChanges.size();iSignChange++)
  {
#ifdef DEBUG_ACTIVITY_DETECTION_ERROR_RATE
	cout << "computing activity detection from " << _signChanges[iSignChange-1] << " to " << _signChanges[iSignChange] << " (" << _signChanges[iSignChange]-_signChanges[iSignChange-1] << ")" << endl;
	cout << "permutation being applied is" << endl;
	Util::print(permutations[_piecesBestPermuationIndexes[iSignChange-1]]);
	cout << endl;
#endif
	res += (_signChanges[iSignChange]-_signChanges[iSignChange-1])*
			computeSERwithoutSolvingAmbiguity(
			  sourceSymbols.block(0,_signChanges[iSignChange-1],N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
			  Util::applyPermutationOnRows(detectedSymbols,permutations[_piecesBestPermuationIndexes[iSignChange-1]],vector<int>(N,1)).block(0,_signChanges[iSignChange-1],N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
			  Util::block(mask,0,_signChanges[iSignChange-1],N,_signChanges[iSignChange]-_signChanges[iSignChange-1]));
  }

  res /= sourceSymbols.cols();
  
  return res;
}

void CDMASystem::BeforeEndingAlgorithm()
{
  SMCSystem::BeforeEndingAlgorithm();

  double peActivity = computeActivityDetectionErrorRate(symbols.block(0,preambleLength,N,frameLength),detectedSymbols);
  
  // the activity detection error probability is accumulated
  _presentFramePeActivityDetection(iSNR,iAlgorithm) = peActivity;
}

void CDMASystem::OnlyOnce()
{
	SMCSystem::OnlyOnce();

	_presentFramePeActivityDetection = MatrixXd::Zero(SNRs.size(),algorithms.size());
}

bool CDMASystem::isChannelOk(const MIMOChannel * const channel)
{
  double thisChannelMatrixMaximumRatio;
  int iChannelMatrix;
  
  // we check if the channel is really bad (severe near-far issues)...
  _maximumRatio = 20*log10(Util::maxCoefficientsRatio(channel->at(preambleLength)));
  
  //...or if any of its coefficients changes sign
  MatrixXd firstSignsMatrix = Util::sign(channel->at(preambleLength));
  
  // all the channel matrices contained in this channel are checked
  for(iChannelMatrix=preambleLength+1;iChannelMatrix<channel->length();iChannelMatrix++)
  {
	// check if how large is the difference of power between coefficients
	thisChannelMatrixMaximumRatio = 20*log10(Util::maxCoefficientsRatio(channel->at(iChannelMatrix)));
	if(thisChannelMatrixMaximumRatio < _maximumRatio)
	  _maximumRatio = thisChannelMatrixMaximumRatio;
	
// 	// check if any coefficient changes sign
	if(Util::sign(channel->at(iChannelMatrix))!=firstSignsMatrix)
	{
// 	  cout << COLOR_PINK << "Coefficients change sign...channel is NOT ok!!" << COLOR_NORMAL << endl;
// 	  return false;
	}
  }
  
  cout << COLOR_WHITE << "the max difference among coefficients in dBs: " << COLOR_NORMAL << _maximumRatio << endl;

#ifdef DEBUG
  cout << channel->nOutputs() << " " << channel->nInputs() << endl;
#endif
  
  if(_maximumRatio > _maximumRatioThresholdInDBs)
  {
	cout << COLOR_PINK << "the max difference among coefficients in dBs is " << COLOR_NORMAL << _maximumRatio << COLOR_PINK << "...channel is NOT ok!!" << COLOR_NORMAL << endl;
	return false;
  }
  
  return true;
}

double CDMASystem::computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns)
{
  if(symbolsDetectionWindowStart!=0)
	throw RuntimeException("CDMASystem::computeSER: this is only implemented when the starting point for SER computing is zero (the beginning of the frame).");
  
  if(detectedSymbols.rows()==0)
  {
	_piecesBestPermuationIndexes = std::vector<uint>(_signChanges.size()-1,0);
	_piecesBestPermutationSigns = std::vector<std::vector<int> >(_signChanges.size()-1,std::vector<int>(N,1));
	return -1.0;
  }

  _piecesBestPermuationIndexes.clear();
  _piecesBestPermutationSigns.clear();  

  double res = 0.0;
  
//   _signChanges = channel->getInputsZeroCrossings(preambleLength+symbolsDetectionWindowStart,frameLength-symbolsDetectionWindowStart);
  //								^
  //								|
  // even though ultimately only symbol vectors from "symbolsDetectionWindowStart" will be taken into account for detection, we don't
  // have to worry about it here, since the mask will take care of that. That being so, "_signChanges" actually contains all the sign changes
  // that occur within the frame
  _signChanges = channel->getInputsZeroCrossings(preambleLength,frameLength);

  uint overallNumberAccountedSymbols = 0;
  
  for(uint iSignChange=1;iSignChange<_signChanges.size();iSignChange++)
  {
	  // we need to find out how many symbols are gonna be taken into account within this subframe
	  uint thisSubframeNumberAccountedSymbols = 0;
	  for(int i=0;i<N;i++)
		for(uint j=_signChanges[iSignChange-1];j<_signChanges[iSignChange];j++)
		  thisSubframeNumberAccountedSymbols += mask[i][j];

	  overallNumberAccountedSymbols += thisSubframeNumberAccountedSymbols;
	  
#ifdef DEBUG_SER
	  cout << "thisSubframeNumberAccountedSymbols = " << thisSubframeNumberAccountedSymbols << endl;
#endif

// 	  res += (_signChanges[iSignChange]-_signChanges[iSignChange-1])*
	  res += thisSubframeNumberAccountedSymbols*
				BaseSystem::computeSER(sourceSymbols.block(0,_signChanges[iSignChange-1],N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
				detectedSymbols.block(0,_signChanges[iSignChange-1],N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
				Util::block(mask,0,_signChanges[iSignChange-1],N,_signChanges[iSignChange]-_signChanges[iSignChange-1]),
				iBestPermutation,bestPermutationSigns);
				
	  // we need to store which the best permutations were along with their corresponding signs
	  _piecesBestPermuationIndexes.push_back(iBestPermutation);
	  _piecesBestPermutationSigns.push_back(bestPermutationSigns);
  }

#ifdef DEBUG_SER
//   cout << "overallNumberAccountedSymbols = " << overallNumberAccountedSymbols << endl;
//   cout << "_piecesBestPermuationIndexes.size() = " << _piecesBestPermuationIndexes.size() << endl;
#endif

//   res /= sourceSymbols.cols();
//   res /= sourceSymbols.cols()*sourceSymbols.rows();
  res /= overallNumberAccountedSymbols;

  return res;
}

double CDMASystem::computeMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices) const
{
  double res = 0.0;

  if(estimatedChannelMatrices.size()==0)
	return -1.0;
  
  if(realChannelMatrices.size()!=estimatedChannelMatrices.size())
	throw RuntimeException("CDMASystem::computeMSE: different number of real channel matrices than estimated channel matrices.");
  
#ifdef DEBUG_MSE
	cout << "realChannelMatrices.size() = " << realChannelMatrices.size() << " estimatedChannelMatrices.size() = " << estimatedChannelMatrices.size() << endl;
	for(int j=0;j<10;j++)
	{
	  cout << "realChannelMatrices[j]" << endl << realChannelMatrices[j] << endl;
	  cout << "estimatedChannelMatrices[j]" << endl << estimatedChannelMatrices[j] << endl;
	}
#endif

//   __signChanges = channel->getInputsZeroCrossings(preambleLength+symbolsDetectionWindowStart,frameLength-symbolsDetectionWindowStart);
  
  // we know that the received channel matrices go from time instant preambleLength+MSEwindowStart (both known parameters of the system) 
  // to the end of the frame. The first matrix of realChannelMatrices/estimatedChannelMatrices corresponds to the time instant:
  uint iChannelMatricesStart = preambleLength+MSEwindowStart;
  
  // we find which interval "iChannelMatricesStart" belongs to
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
  
  res += (_signChanges[iSignChange]-iChannelMatricesStart)*BaseSystem::computeMSE(toCheckRealChannelMatrices,toCheckEstimatedChannelMatrices,
								permutations[_piecesBestPermuationIndexes[iSignChange-1]],_piecesBestPermutationSigns[iSignChange-1]);

  iSignChange++;
  
  for(;iSignChange<_signChanges.size();iSignChange++)
  {
	toCheckRealChannelMatrices = std::vector<MatrixXd>(realChannelMatrices.begin()+_signChanges[iSignChange-1]-iChannelMatricesStart,realChannelMatrices.begin()+_signChanges[iSignChange]-iChannelMatricesStart);
	toCheckEstimatedChannelMatrices = std::vector<MatrixXd>(estimatedChannelMatrices.begin()+_signChanges[iSignChange-1]-iChannelMatricesStart,estimatedChannelMatrices.begin()+_signChanges[iSignChange]-iChannelMatricesStart);
	
	res += (_signChanges[iSignChange]-iChannelMatricesStart)*BaseSystem::computeMSE(toCheckRealChannelMatrices,toCheckEstimatedChannelMatrices,
								permutations[_piecesBestPermuationIndexes[iSignChange-1]],_piecesBestPermutationSigns[iSignChange-1]);

#ifdef DEBUG_MSE
	cout << "computing MSE between " << _signChanges[iSignChange-1] << " and " << _signChanges[iSignChange] << " ( " << (_signChanges[iSignChange]-_signChanges[iSignChange-1]) << " matrices)" << endl;
#endif
  }

  res /= realChannelMatrices.size();

  return res;
}