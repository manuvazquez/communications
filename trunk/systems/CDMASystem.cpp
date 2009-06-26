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

#define DEBUG

CDMASystem::CDMASystem(): BaseSystem(),ARcoefficients(1)
{
    nParticles = 192;
    resamplingRatio = 0.9;
    
    if(m!=1)
        throw RuntimeException("CDMASystem::CDMASystem: channel is not flat.");
    
    // spreading spreadingCodes for the users are generated randomly
    _spreadingCodes = StatUtil::RandnMatrix(L,N,0.0,1.0);
    _spreadingCodes = Util::sign(_spreadingCodes);
    
#ifdef DEBUG
    cout << "generated spreadingCodes..." << endl << _spreadingCodes;
#endif
    
    // AR process parameters
    ARcoefficients[0] = 0.99999;
    ARvariance=0.0001;    
    
    // a flat power profile is generated. Notice:
    //      i) that m should be 1, otherwise an exception would have been thrown
    //     ii) we only need to generate a coefficient per user, i.e., a Nx1 vector
    powerProfile = new FlatPowerProfile(N,1,m,1.0);

    // particle filtering
    ResamplingCriterion criterioRemuestreo(resamplingRatio);
    algoritmoRemuestreo = new ResidualResamplingAlgorithm(criterioRemuestreo);

    userPersistenceProb = 0.8;
    userActivityProb = 0.2;
    userPriorProb = 0.5;
//     exit(0);
}


CDMASystem::~CDMASystem()
{
    delete powerProfile;
    delete algoritmoRemuestreo;
}


void CDMASystem::AddAlgorithms()
{
#ifdef DEBUG
    cout << "observations are" << endl << observations;
#endif
}

void CDMASystem::BeforeEndingAlgorithm(int iAlgorithm)
{
    BaseSystem::BeforeEndingAlgorithm(iAlgorithm);
}

void CDMASystem::BeforeEndingFrame(int iFrame)
{
    BaseSystem::BeforeEndingFrame(iFrame);
}

void CDMASystem::BuildChannel()
{
    channel = new ARMultiuserCDMAchannel(symbols.cols(),_spreadingCodes,ARprocess(powerProfile->GenerateChannelMatrix(randomGenerator),ARcoefficients,ARvariance));
    
    // when users are not transmitting, their symbols are zero
    _usersActivity.resize(symbols.rows(),symbols.cols());
    
    // at the first time instant the prior probability is used to decide which users are active
//     for(uint iUser=0;iUser<symbols.rows();iUser++)
      
  
    for(uint iTime=preambleLength+trainSeqLength+1;iTime<symbols.cols();iTime++)
        for(uint iUser=0;iUser<symbols.rows();iUser++)
        {
            
        }
}

