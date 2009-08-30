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
#ifndef CDMASYSTEM_H
#define CDMASYSTEM_H

#include <SMCSystem.h>
#include <FlatPowerProfile.h>
#include <ARMultiuserCDMAchannel.h>
#include <CDMAKalmanEstimator.h>
#include <CDMAunknownActiveUsersSISopt.h>
#include <CDMAKnownChannelChannelMatrixEstimator.h>
#include <KnownFlatChannelOptimalAlgorithm.h>
#include <KnownFlatChannelAndActiveUsersOptimalAlgorithm.h>
#include <UsersActivityDistribution.h>
#include <UnknownActiveUsersLinearFilterBasedSMCAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/
class CDMASystem : public SMCSystem
{
protected:
    MatrixXd _spreadingCodes;
    
    // _usersActivity(i,j) = 1.0 if the i-th user is active at time j
    vector<vector<bool> > _usersActivity;
    
    double userPersistenceProb,newActiveUserProb,userPriorProb;
    
    CDMAKalmanEstimator *cdmaKalmanEstimator;
    CDMAKnownChannelChannelMatrixEstimator *cdmaKnownChannelChannelMatrixEstimator;
    
    MMSEDetector *mmseDetector;
    
    UsersActivityDistribution usersActivityPdf;    
    
    virtual void AddAlgorithms();
    virtual void BeforeEndingFrame(int iFrame);
    virtual void BuildChannel();    
public:
    CDMASystem();

    ~CDMASystem();

};

#endif
