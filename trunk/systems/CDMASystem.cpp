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

// #define DEBUG

CDMASystem::CDMASystem(): SMCSystem()
{
    if(m!=1)
        throw RuntimeException("CDMASystem::CDMASystem: channel is not flat.");
    
    // spreading spreadingCodes for the users are generated randomly
    _spreadingCodes = StatUtil::RandnMatrix(L,N,0.0,1.0);
    _spreadingCodes = Util::sign(_spreadingCodes);
    
#ifdef DEBUG
    cout << "generated spreadingCodes..." << endl << _spreadingCodes;
#endif
    
    // AR process parameters
    ARcoefficients = vector<double>(2);
    ARcoefficients[0] = 0.59999;
    ARcoefficients[1] = 0.39999;
    ARvariance=0.0001;    
    
    // a flat power profile is generated. Notice:
    //      i) that m should be 1, otherwise an exception would have been thrown
    //     ii) we only need to generate a coefficient per user, i.e., a 1xN vector
    powerProfile = new FlatPowerProfile(1,N,m,1.0);

    userPersistenceProb = 0.8;
    newActiveUserProb = 0.2;
    userPriorProb = 0.5;
    
    userPersistenceProb = 0.8;
    newActiveUserProb = 0.2;
    userPriorProb = 1.0;    
    
    cdmaKalmanEstimator = new CDMAKalmanEstimator(powerProfile->means(),powerProfile->variances(),ARcoefficients,ARvariance,_spreadingCodes);
    cdmaKnownChannelChannelMatrixEstimator = NULL;
}


CDMASystem::~CDMASystem()
{
    delete powerProfile;
    delete cdmaKalmanEstimator;
    delete cdmaKnownChannelChannelMatrixEstimator;
}

void CDMASystem::AddAlgorithms()
{
#ifdef DEBUG
    cout << "observations are" << endl << observations;
#endif
    algorithms.push_back(new CDMAunknownActiveUsersSISopt ("CDMA SIS-opt",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKalmanEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),userPersistenceProb,newActiveUserProb,userPriorProb));
       
//     // the channel is different in each frame, so the estimator that knows the channel must be rebuilt every frame
//     delete cdmaKnownChannelChannelMatrixEstimator;
//     cdmaKnownChannelChannelMatrixEstimator = new CDMAKnownChannelChannelMatrixEstimator(channel,preambleLength,N,_spreadingCodes);
//  
//     algorithms.push_back(new CDMAunknownActiveUsersSISopt ("CDMA SIS-opt (known channel)",*alphabet,L,1,N,iLastSymbolVectorToBeDetected,m,cdmaKnownChannelChannelMatrixEstimator,preamble,d,nParticles,algoritmoRemuestreo,powerProfile->means(),powerProfile->variances(),userPersistenceProb,newActiveUserProb,userPriorProb));
}

void CDMASystem::BeforeEndingAlgorithm(int iAlgorithm)
{
    SMCSystem::BeforeEndingAlgorithm(iAlgorithm);
}

void CDMASystem::BeforeEndingFrame(int iFrame)
{
    SMCSystem::BeforeEndingFrame(iFrame);
}

void CDMASystem::BuildChannel()
{
  
#ifdef DEBUG
  cout << "symbols before" << endl << symbols;
#endif
    
    // when users are not transmitting, their symbols are zero
    _usersActivity = vector<vector<bool> >(symbols.rows(),vector<bool>(frameLength));
    
    tVector userActivePriorPdf(2);
    userActivePriorPdf(0) = 1.0 - userPriorProb;
    userActivePriorPdf(1) = userPriorProb;
    
    tVector userActiveGivenItWasPdf(2);
    userActiveGivenItWasPdf(0) = 1.0 - userPersistenceProb;
    userActiveGivenItWasPdf(1) = userPersistenceProb;    
    
    tVector userActiveGivenItWasNotPdf(2);
    userActiveGivenItWasNotPdf(0) = 1.0 - newActiveUserProb;
    userActiveGivenItWasNotPdf(1) = newActiveUserProb;        
    
    vector<int> usersActive = StatUtil::discrete_rnd(symbols.rows(),userActivePriorPdf);
    
    // at the first time instant the prior probability is used to decide which users are active
    for(uint iUser=0;iUser<static_cast<uint>(symbols.rows());iUser++)
    {
        _usersActivity[iUser][trainSeqLength] = bool(usersActive[iUser]);
        symbols(iUser,preambleLength+trainSeqLength) = double(_usersActivity[iUser][trainSeqLength])*symbols(iUser,preambleLength+trainSeqLength);
        isSymbolAccountedForDetection[iUser][trainSeqLength] = _usersActivity[iUser][trainSeqLength];
    }
      
    // set of active users evolves according to the given probabilities
    for(uint iTime=trainSeqLength+1;iTime<frameLength;iTime++)    
        for(uint iUser=0;iUser<symbols.rows();iUser++)
        {   
            // the user was active in the last time instant
            if(_usersActivity[iUser][iTime-1]) 
                _usersActivity[iUser][iTime]= bool(StatUtil::discrete_rnd(userActiveGivenItWasPdf));
            // the user was NOT active in the last time instant
            else
                _usersActivity[iUser][iTime]= bool(StatUtil::discrete_rnd(userActiveGivenItWasNotPdf));
            
            symbols(iUser,preambleLength+iTime) = symbols(iUser,preambleLength+iTime)*double(_usersActivity[iUser][iTime]);
            isSymbolAccountedForDetection[iUser][iTime] = _usersActivity[iUser][iTime];
        }
            
#ifdef DEBUG
    Util::print(usersActive);
    cout << "users activity at time 0" << endl;
    Util::print(_usersActivity);
#endif            
                
#ifdef DEBUG
    cout << "symbols after" << endl << symbols;
#endif    
    
    channel = new ARMultiuserCDMAchannel(symbols.cols(),_spreadingCodes,ARprocess(powerProfile->generateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
}
