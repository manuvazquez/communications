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
#ifndef UNKNOWNACTIVEUSERSLINEARFILTERBASEDSMCALGORITHM_H
#define UNKNOWNACTIVEUSERSLINEARFILTERBASEDSMCALGORITHM_H

#include <SMCAlgorithm.h>

/*!
An SMC algorithm based on linear filters that for a system whose users are not permanently active (at the moment only for flat channels)

    @author Manu <manu@rustneversleeps>
*/

#include <ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers.h>
#include <KalmanEstimator.h>
#include <UsersActivityDistribution.h>

class UnknownActiveUsersLinearFilterBasedSMCAlgorithm : public SMCAlgorithm
{
public:
    UnknownActiveUsersLinearFilterBasedSMCAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, LinearDetector *linearDetector, MatrixXd preamble, uint smoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm, const MatrixXd& channelMatrixMean, const MatrixXd& channelMatrixVariances, const std::vector<UsersActivityDistribution> usersActivityPdfs);

    ~UnknownActiveUsersLinearFilterBasedSMCAlgorithm();

protected:
    LinearDetector *_linearDetector;
    const std::vector<UsersActivityDistribution> _usersActivityPdfs; /// objects describing the pdf of the users activity    

    virtual void initializeParticles();
    virtual void process(const MatrixXd& observations, vector< double > noiseVariances);

private:
	double probSymbol(tSymbol symbol, UsersActivityDistribution userActivityDistribution) const;
	double probSymbolGivenPreviousActivity(tSymbol symbol,bool previousActivity,UsersActivityDistribution userActivityDistribution) const;
};

#endif
