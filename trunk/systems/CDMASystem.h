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

#include <BaseSystem.h>
#include <FlatPowerProfile.h>
#include <ARMultiuserCDMAchannel.h>
#include <ResamplingAlgorithm.h>
#include <ResidualResamplingAlgorithm.h>
// #include <CDMAKalmanEstimator.h>

/**
    @author Manu <manu@rustneversleeps>
*/
class CDMASystem : public BaseSystem
{
public:
    CDMASystem();

    ~CDMASystem();

protected:
    tMatrix _spreadingCodes;
    
    // _usersActivity(i,j) = 1.0 if the i-th user is active at time j
    vector<vector<bool> > _usersActivity;
    
    double userPersistenceProb,newActiveUserProb,userPriorProb;
    
    std::vector<double> ARcoefficients;
    double ARvariance;

    int nParticles;
    double resamplingRatio;
    ResamplingAlgorithm *algoritmoRemuestreo;
    
//     CDMAKalmanEstimator *cdmaKalmanEstimator;

    virtual void AddAlgorithms();
    virtual void BeforeEndingAlgorithm(int iAlgorithm);
    virtual void BeforeEndingFrame(int iFrame);
    virtual void BuildChannel();

};

#endif
