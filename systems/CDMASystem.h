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
#include <KnownSymbolsKalmanBasedChannelEstimatorAlgorithm.h>
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
    
    double _userPersistenceProb,_newActiveUserProb,_userPriorProb;
	
    double _velocity,_carrierFrequency,_symbolRate,_T;
    
    CDMAKalmanEstimator *_cdmaKalmanEstimator;
    CDMAKnownChannelChannelMatrixEstimator *_cdmaKnownChannelChannelMatrixEstimator;
	
    MMSEDetector *_mmseDetector;
    
	std::vector<UsersActivityDistribution> _usersActivityPdfs;
	
	// it stores the maximum ratio among the coefficients of a channel matrix for every frame
	vector<double> _maxCoefficientsRatiosInDBs;
	
	// this is used in two different methods (though only computed in one of them...)
	double _maximumRatio;
	
	double _maximumRatioThresholdInDBs;
	
	MatrixXd _presentFramePeActivityDetection;
	vector<MatrixXd> _peActivityDetectionFrames;
	
	int _nSurvivors;

    bool _adjustSurvivorsFromParticlesNumber;
    bool _adjustParticlesNumberFromSurvivors;
	
	// needed for solving ambiguity when computing MSE
	
	// each element of the vector (a vector itself) stores the signs of the best permutation corresponding
	// to one of the pieces into which the frame was dividing in order to compute the SER
	std::vector<std::vector<int> > _piecesBestPermutationSigns;
	
	// each element of the vector stores the index of the best permutation corresponding to one of the
	// pieces into which the frame was dividing in order to compute the SER
	std::vector<uint> _piecesBestPermuationIndexes;
	
	// a vector that stores the indexes within the frame (from 0 to "frameLength") where a sign change takes place
	// it ALWAYS includes the first index (0) and the last ("frameLength")
	std::vector<uint> _signChanges;
	
    virtual void addAlgorithms();
	virtual void beforeEndingAlgorithm();
    virtual void beforeEndingFrame();
    virtual void buildChannel();
	virtual void onlyOnce();
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
	double computeActivityDetectionErrorRate(MatrixXd sourceSymbols,MatrixXd detectedSymbols) const;
	
	virtual double computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns);
	virtual double computeMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices) const;

  private:
	  
	//! It checks the channel for high differences in the coefficients power (near-far problem)
	/*!
	  \param channel a pointer to the channel to check
	  \return true if maximum difference among coefficients (in dBs) is higher than the threshold \ref _maximumRatioThresholdInDBs
	*/
	bool isChannelOk(const MIMOChannel * const channel);
};

#endif
