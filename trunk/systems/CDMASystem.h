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
#include <CDMAKalmanEstimator.h>
#include <CDMAKnownChannelChannelMatrixEstimator.h>
#include <UsersActivityDistribution.h>

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
	
    double velocity,carrierFrequency,symbolRate,T;
    
    CDMAKalmanEstimator *cdmaKalmanEstimator;
    CDMAKnownChannelChannelMatrixEstimator *cdmaKnownChannelChannelMatrixEstimator;
    
    MMSEDetector *mmseDetector;
    
	std::vector<UsersActivityDistribution> usersActivityPdfs;
	
	// it stores the maximum ratio among the coefficients of a channel matrix for every frame
	vector<double> _maxCoefficientsRatiosInDBs;
	
	// this is used in two different methods (though only computed in one of them...)
	double _maximumRatio;
	
	double maximumRatioThresholdInDBs;
	
	MatrixXd presentFramePeActivityDetection;
	vector<MatrixXd> peActivityDetectionFrames;
	
	int nSurvivors;

    bool adjustSurvivorsFromParticlesNumber;
    bool adjustParticlesNumberFromSurvivors;
	
    virtual void AddAlgorithms();
	virtual void BeforeEndingAlgorithm(int iAlgorithm);
    virtual void BeforeEndingFrame(int iFrame);
    virtual void BuildChannel();
	virtual void OnlyOnce();
public:
    CDMASystem();

    ~CDMASystem();

	bool areSequencesOrthogonal(const MatrixXd &spreadingCodes);
	
	/*!
	  It computes the probability of detecting that a user is transmitting when it's not or the other way around. It relies in \ref computeSER
	  \param sourceSymbols are the actual transmitted symbols
	  \param detectedSymbols are the symbols detected
	  \return the probability of misdetection the activity of a user
	*/
	double computeActivityDetectionER(MatrixXd sourceSymbols,MatrixXd detectedSymbols);
	static bool isUserActive(const tSymbol symbol) { return symbol!=0.0;}
  private:
	bool isChannelOk(const MIMOChannel * const channel);
};

#endif
