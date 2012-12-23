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
#include <CDMARLSEstimator.h>
#include <KnownSymbolsKalmanBasedChannelEstimatorAlgorithm.h>
#include <UsersActivityDistribution.h>

// #define ESTIMATE_CHANNEL_TRANSITION_PROBABILITIES

// when dealing with a high number of frames, the users' activiy and spreading codes for every frame may amount to a lot of space in memory/disk
// #define KEEP_EVERY_FRAME_USER_ACTIVITY
// #define KEEP_EVERY_FRAME_SPREADING_CODES

/**
	@author Manu <manu@rustneversleeps>
*/
class CDMASystem : public SMCSystem
{
	
#ifdef ESTIMATE_CHANNEL_TRANSITION_PROBABILITIES
private:
	uint channelCoeffToCell(double coeff, const std::vector<double> &grid, double gridStep) const;
#endif
	
protected:
    MatrixXd _spreadingCodes;
    
    // _usersActivity(i,j) = 1.0 if the i-th user is active at time j
    vector<vector<bool> > _usersActivity;
	
#ifdef KEEP_EVERY_FRAME_USER_ACTIVITY
	std::vector<std::vector<std::vector<bool> > > _everyFrameUsersActivity;
#endif
    
    double _userPersistenceProb,_newActiveUserProb,_userPriorProb;
	
    CDMAKalmanEstimator *_cdmaKalmanEstimator;
	CDMAKnownChannelChannelMatrixEstimator *_cdmaKnownChannelChannelMatrixEstimator;
	CDMARLSEstimator *_cdmaRLSEstimator;
	
    MMSEDetector *_mmseDetector;
    
	std::vector<UsersActivityDistribution> _usersActivityPdfs;
	
	MatrixXd _presentFramePeActivityDetection;
	vector<MatrixXd> _peActivityDetectionFrames;
	
	uint _nSurvivors;

    bool _adjustSurvivorsFromParticlesNumber;
    bool _adjustParticlesNumberFromSurvivors;
	
	// needed for solving ambiguity when computing MSE
	
	// each element of the vector (a vector itself) stores the signs of the best permutation corresponding to one of the pieces into which the frame was dividing in order to compute the SER
	std::vector<std::vector<int> > _piecesBestPermutationSigns;
	
	// each element of the vector stores the index of the best permutation corresponding to one of the pieces into which the frame was dividing in order to compute the SER
	std::vector<uint> _piecesBestPermuationIndexes;
	
	// a vector that stores the indexes within the frame (from 0 to "frameLength") where a sign change takes place. It ALWAYS includes the first index (0) and the last ("frameLength")
	std::vector<uint> _signChanges;
	
	/**
	 * @brief the user of interest, for which is computed the BER and whose SNR is fixed.
	 **/
	uint _iUserOfInterest;
	
	/**
	 * @brief if the signal to interference ratio is below this threshold the channel is discarded
	 **/
	int _minSignalToInterferenceRatio;
	
#ifdef KEEP_EVERY_FRAME_SPREADING_CODES
	std::vector<MatrixXd> _everyFrameSpreadingCodes;
#endif
	
	/**
	 * @brief number of sign changes that occur in the channel estimate of all the algorithms for every SNR and every frame
	 **/
	std::vector<std::vector<std::vector<uint> > > _everyFrameNumberSignChanges;
	
	std::vector<std::vector<uint> > _thisFrameNumberSignChanges;
	
	std::string _maskUsedToComputeTheSER;
	
	double _forgettingFactor;
	
// 	std::vector<double> _grid,_grid20,_grid30,_grid50;
// 	MatrixXd _channelTransitionProbabilities,_channelTransitionProbabilities20,_channelTransitionProbabilities30,_channelTransitionProbabilities50;
	
	std::vector<std::vector<double> > _grids;
	std::vector<std::string> _channelTransitionProbabilitiesFileNames;
	std::vector<std::string> _channelMarginalProbabilitiesFileNames;
	std::vector<uint> _nCells;
	std::vector<MatrixXd> _estimatedChannelTransitionProbabilities;
	std::vector<MatrixXd> _estimatedChannelMarginalProbabilities;
	
#ifdef ESTIMATE_CHANNEL_TRANSITION_PROBABILITIES
	std::vector<double> _gridSteps;
#endif
	
// 	std::string _channelTransitionProbabilitiesFileName;
	
// #ifdef ESTIMATE_CHANNEL_TRANSITION_PROBABILITIES
// 	MatrixXd _estimatedChannelTransitionProbabilities;
// 	VectorXd _estimatedChannelMarginalProbabilities;
// #endif
	
    virtual void addAlgorithms();
	virtual void beforeEndingAlgorithm();
	virtual void storeFrameResults();
	virtual void saveFrameResults();
    virtual void buildSystemSpecificVariables();
	virtual void onlyOnce();
	
	/**
	 * @brief it initializes the variables that account for the splitting of the frame so that: i) the frame is not split ii) the best permutation is the first and its associated sign is +1
	 *
	 * @return void
	 **/
	void resetFramePieces();
	
#ifdef ESTIMATE_CHANNEL_TRANSITION_PROBABILITIES
	void accountForEstimatedChannelTransitionProbabilities(const MIMOChannel * const channel);
#endif
	
	void readChannelCoefficientsGridsParametersFromXML(xml_node<> *parentNode,std::string xmlName);
	
public:
    CDMASystem();

    ~CDMASystem();

	/*!
	  It computes the probability of detecting that a user is transmitting when it's not or the other way around. It relies in \ref computeSER
	  \param sourceSymbols are the actual transmitted symbols
	  \param detectedSymbols are the symbols detected
	  \return the probability of misdetection the activity of a user
	*/
	double computeActivityDetectionErrorRate(MatrixXd sourceSymbols,MatrixXd detectedSymbols) const;
	
	double computeSelectedUsersActivityDetectionErrorRate(MatrixXd sourceSymbols,MatrixXd detectedSymbols) const;
	
	virtual double computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns);
	virtual double computeSelectedUsersSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns);
	virtual double computeMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices,const std::vector<bool> &mask) const;
	virtual double computeSelectedUsersMSE(const vector<MatrixXd> &realChannelMatrices,const vector<MatrixXd> &estimatedChannelMatrices,const std::vector<bool> &mask) const;
	
	virtual Noise *buildNoise() const;
	virtual MIMOChannel *buildChannel();
};

#endif
