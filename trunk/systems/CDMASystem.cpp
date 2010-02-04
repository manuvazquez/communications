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

#include <math.h>

// #define PRINT_INFO

CDMASystem::CDMASystem(): SMCSystem()
,userPersistenceProb(0.8),newActiveUserProb(0.2),userPriorProb(1.0)
// ,userPersistenceProb(1.0),newActiveUserProb(0.2),userPriorProb(1.0)
,usersActivityPdf(userPersistenceProb,newActiveUserProb,userPriorProb)
{
    if(m!=1)
        throw RuntimeException("CDMASystem::CDMASystem: channel is not flat.");
    
    // spreading spreadingCodes for the users are generated randomly
    _spreadingCodes = StatUtil::randnMatrix(L,N,0.0,1.0);
    _spreadingCodes = Util::sign(_spreadingCodes);
	
	// the spreading codes are normalized
// 	_spreadingCodes /= sqrt(L);
    
#ifdef PRINT_INFO
    cout << "generated spreadingCodes..." << endl << _spreadingCodes;
#endif
    
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
    algorithms.push_back(new KnownFlatChannelOptimalAlgorithm ("CDMA optimal (known channel)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,*channel,preambleLength));
    
    algorithms.push_back(new KnownFlatChannelAndActiveUsersOptimalAlgorithm ("CDMA optimal (known channel & active users)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,*channel,preambleLength,_usersActivity));    
       
    // the channel is different in each frame, so the estimator that knows the channel must be rebuilt every frame
    delete cdmaKnownChannelChannelMatrixEstimator;
    cdmaKnownChannelChannelMatrixEstimator = new CDMAKnownChannelChannelMatrixEstimator(channel,preambleLength,N,_spreadingCodes);
 
//     algorithms.push_back(new CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge ("CDMA SIS-opt with no knowledge of users activity pdf (known channel)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKnownChannelChannelMatrixEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means_eigen(),powerProfile->variances_eigen()));
    
    algorithms.push_back(new CDMAunknownActiveUsersSISoptWithNoUsersActivityKnowledge ("CDMA SIS-opt with no knowledge of users activity pdf",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances()));

    algorithms.push_back(new CDMAunknownActiveUsersSISopt ("CDMA SIS-opt",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),usersActivityPdf));	
	
    algorithms.push_back(new UnknownActiveUsersLinearFilterBasedSMCAlgorithm ("CDMA SIS Linear Filters",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,mmseDetector,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),usersActivityPdf));            
}

void CDMASystem::BeforeEndingFrame(int iFrame)
{
    SMCSystem::BeforeEndingFrame(iFrame);
    Util::matrixToOctaveFileStream(_spreadingCodes,"spreadingCodes",f);
}

void CDMASystem::BuildChannel()
{
    // when users are not transmitting, their symbols are zero
    _usersActivity = vector<vector<bool> >(symbols.rows(),vector<bool>(frameLength));
    
    // at the first time instant the prior probability is used to decide which users are active
    for(uint iUser=0;iUser<static_cast<uint>(symbols.rows());iUser++)
    {
        _usersActivity[iUser][trainSeqLength] = usersActivityPdf.sampleFromPrior();        
        symbols(iUser,preambleLength+trainSeqLength) = double(_usersActivity[iUser][trainSeqLength])*symbols(iUser,preambleLength+trainSeqLength);
        isSymbolAccountedForDetection[iUser][trainSeqLength] = _usersActivity[iUser][trainSeqLength];
    }
      
    // set of active users evolves according to the given probabilities
    for(int iTime=trainSeqLength+1;iTime<frameLength;iTime++)    
        for(int iUser=0;iUser<symbols.rows();iUser++)
        {   
            _usersActivity[iUser][iTime] = usersActivityPdf.sampleGivenItWas(_usersActivity[iUser][iTime-1]);             
            symbols(iUser,preambleLength+iTime) = symbols(iUser,preambleLength+iTime)*double(_usersActivity[iUser][iTime]);
            isSymbolAccountedForDetection[iUser][iTime] = _usersActivity[iUser][iTime];
        }
            
#ifdef PRINT_INFO
    cout << "symbols after generating users activity" << endl << symbols;
#endif    
    
    channel = new ARMultiuserCDMAchannel(symbols.cols(),_spreadingCodes,ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
}
