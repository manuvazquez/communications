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

#define PRINT_CODES_INFO
// #define PRINT_ACTIVITY_SAMPLING

// #define PRINT_INFO

CDMASystem::CDMASystem(): SMCSystem()
// ,userPersistenceProb(0.99),newActiveUserProb(0.01),userPriorProb(0.5)
,userPersistenceProb(0.8),newActiveUserProb(0.2),userPriorProb(1.0)
// ,userPersistenceProb(1.0),newActiveUserProb(0.2),userPriorProb(1.0)
,usersActivityPdfs(N,UsersActivityDistribution(userPersistenceProb,newActiveUserProb,userPriorProb))
// ,maximumRatioThresholdInDBs(15)
,maximumRatioThresholdInDBs(20)
{
    if(m!=1)
        throw RuntimeException("CDMASystem::CDMASystem: channel is not flat.");

	// first users starts transmitting something
	usersActivityPdfs[0].setApriori(1.0);
	
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

	nSurvivors = 2;
// 	nSurvivors = 8;
// 	nSurvivors = 20;

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
    
    cdmaKalmanEstimator = new CDMAKalmanEstimator(powerProfile->means(),powerProfile->variances(),ARcoefficients,ARvariance,_spreadingCodes);
    cdmaKnownChannelChannelMatrixEstimator = NULL;
    
    mmseDetector = new MMSEDetector(L,N,alphabet->variance(),N);

	// bessel channel parameters
    velocity = 50/3.6; // (m/s)
    carrierFrequency = 2e9; // (Hz)
    symbolRate = 500e3; // (Hz)

    T = 1.0/symbolRate; // (s)
	
	_maxCoefficientsRatiosInDBs.reserve(nFrames);
	
	peActivityDetectionFrames.reserve(nFrames);
	
	// adjusting the number of particles from that of the survivors or the other way around
	adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;
	
    // check the adjustments for particle and survivor numbers
    if(adjustParticlesNumberFromSurvivors && adjustSurvivorsFromParticlesNumber)
        throw RuntimeException("CDMASystem::CDMASystem: \"adjustParticlesNumberFromSurvivors\" and \"adjustSurvivorsFromParticlesNumber\" shouldn't be true at the same time.");

    if(adjustParticlesNumberFromSurvivors)
    {
		cout << "Number of particles adjusted from " << nParticles;
        nParticles = (int)pow((double)alphabet->length()+1,N)*nSurvivors;
        cout << " to " << nParticles << endl;
    }

    if(adjustSurvivorsFromParticlesNumber)
    {
        cout << "Number of survivors adjusted from " << nSurvivors;
        nSurvivors = int(ceil(double(nParticles)/pow((double)alphabet->length()+1,double(N))));
        cout << " to " << nSurvivors << endl;
    }
}


CDMASystem::~CDMASystem()
{
    delete powerProfile;
    delete cdmaKalmanEstimator;
    delete cdmaKnownChannelChannelMatrixEstimator;
    delete mmseDetector;
}

void CDMASystem::AddAlgorithms()
{
    algorithms.push_back(new KnownFlatChannelOptimalAlgorithm ("CDMA optimal with known channel BUT no knowledge of users activity probabilities)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,*channel,preambleLength));
    
    algorithms.push_back(new KnownFlatChannelAndActiveUsersOptimalAlgorithm ("CDMA optimal (known channel and active users)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,*channel,preambleLength,_usersActivity));    
       
    // the channel is different in each frame, so the estimator that knows the channel must be rebuilt every frame
    delete cdmaKnownChannelChannelMatrixEstimator;
    cdmaKnownChannelChannelMatrixEstimator = new CDMAKnownChannelChannelMatrixEstimator(channel,preambleLength,N,_spreadingCodes);
     
//     algorithms.push_back(new CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge ("CDMA SIS-opt with no knowledge of users activity pdf",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances()));

//     algorithms.push_back(new CDMAunknownActiveUsersSISopt ("CDMA SIS-opt",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),usersActivityPdfs));	
	
//     algorithms.push_back(new UnknownActiveUsersLinearFilterBasedSMCAlgorithm ("CDMA SIS Linear Filters",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,mmseDetector,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),usersActivityPdfs));
	
	algorithms.push_back(new ViterbiAlgorithmWithAprioriProbabilities("Viterbi (known channel)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d,usersActivityPdfs));
	
	algorithms.push_back(new PSPAlgorithmWithAprioriProbabilities("PSP",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,preamble,d,iLastSymbolVectorToBeDetected+d,nSurvivors,usersActivityPdfs));
}

void CDMASystem::BeforeEndingFrame(int iFrame)
{
    SMCSystem::BeforeEndingFrame(iFrame);
	
	// the maximum ratio of this frame is added to the vector
	_maxCoefficientsRatiosInDBs.push_back(_maximumRatio);
	
    peActivityDetectionFrames.push_back(presentFramePeActivityDetection);
    Util::matricesVectorToOctaveFileStream(peActivityDetectionFrames,"peActivityDetectionFrames",f);
	
	Util::scalarsVectorToOctaveFileStream(_maxCoefficientsRatiosInDBs,"maxCoefficientsRatiosInDBs",f);
    Util::matrixToOctaveFileStream(_spreadingCodes,"spreadingCodes",f);
	Util::scalarToOctaveFileStream(maximumRatioThresholdInDBs,"maximumRatioThresholdInDBs",f);
	Util::scalarToOctaveFileStream(nSurvivors,"nSurvivors",f);
	Util::scalarToOctaveFileStream(userPersistenceProb,"userPersistenceProb",f);
	Util::scalarToOctaveFileStream(newActiveUserProb,"newActiveUserProb",f);
	Util::scalarToOctaveFileStream(userPriorProb,"userPriorProb",f);
}

void CDMASystem::BuildChannel()
{
    // when users are not transmitting, their symbols are zero
    _usersActivity = vector<vector<bool> >(symbols.rows(),vector<bool>(frameLength));
    
    // at the first time instant the prior probability is used to decide which users are active
    for(uint iUser=0;iUser<static_cast<uint>(symbols.rows());iUser++)
    {
        _usersActivity[iUser][trainSeqLength] = usersActivityPdfs[iUser].sampleFromPrior();
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
            _usersActivity[iUser][iTime] = usersActivityPdfs[iUser].sampleGivenItWas(_usersActivity[iUser][iTime-1]);             
            symbols(iUser,preambleLength+iTime) = symbols(iUser,preambleLength+iTime)*double(_usersActivity[iUser][iTime]);
            isSymbolAccountedForDetection[iUser][iTime] = _usersActivity[iUser][iTime];
        }
            
#ifdef PRINT_INFO
    cout << "symbols after generating users activity" << endl << symbols << endl;
#endif    

	double thisChannelMatrixMaximumRatio;

	do
	{
	  delete channel;

// 	  channel = new MultiuserCDMAchannel(new ARchannel(N,1,m,symbols.cols(),ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance)),_spreadingCodes);
	  
// 	  channel = new MultiuserCDMAchannel(new TimeInvariantChannel(powerProfile->nInputs(),powerProfile->nOutputs(),m,symbols.cols(),MatrixXd::Ones(powerProfile->nOutputs(),powerProfile->nInputs())),_spreadingCodes);
	  
	  channel = new MultiuserCDMAchannel(new BesselChannel(N,1,m,symbols.cols(),velocity,carrierFrequency,T,*powerProfile),_spreadingCodes);	  

	  // we check if the channel is really bad (severe near-far issues)
	  _maximumRatio = 20*log10(Util::maxCoefficientsRatio(channel->at(preambleLength)));
	  
	  // all the channel matrices contained in this channel are checked
	  for(int i=preambleLength+1;i<channel->length();i++)
	  {
		thisChannelMatrixMaximumRatio = 20*log10(Util::maxCoefficientsRatio(channel->at(i)));
		if(thisChannelMatrixMaximumRatio < _maximumRatio)
		  _maximumRatio = thisChannelMatrixMaximumRatio;
	  }
	  
	  cout << "the max difference among coefficients in dBs: " << _maximumRatio << endl;
	
	} while(_maximumRatio>maximumRatioThresholdInDBs);
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

double CDMASystem::computeActivityDetectionER(MatrixXd sourceSymbols, MatrixXd detectedSymbols)
{
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
	
    for(int iTime=0;iTime<symbolsDetectionWindowStart;iTime++)
        for(int iInput=0;iInput<N;iInput++)
            isSymbolAccountedForDetection[iInput][iTime] = false;        
  

	// in order to compute the probability of activity detection it makes no difference the symbol detected
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

	  return computeSER(sourceSymbols,detectedSymbols,mask);
}

void CDMASystem::BeforeEndingAlgorithm(int iAlgorithm)
{
	SMCSystem::BeforeEndingAlgorithm(iAlgorithm);

	double peActivity = computeActivityDetectionER(symbols.block(0,preambleLength,N,frameLength),detectedSymbols);
	
    // the activity detection error probability is accumulated
    presentFramePeActivityDetection(iSNR,iAlgorithm) = peActivity;	
}

void CDMASystem::OnlyOnce()
{
	SMCSystem::OnlyOnce();

	presentFramePeActivityDetection = MatrixXd::Zero(SNRs.size(),algorithms.size());
}

bool CDMASystem::isChannelOk(const ChannelMatrixEstimator * const channel)
{
  return true;
}