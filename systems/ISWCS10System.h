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
#ifndef ISWCS10SYSTEM_H
#define ISWCS10SYSTEM_H

#include <ChannelOrderEstimationSystem.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <RLSEstimator.h>
#include <RMMSEDetector.h>
#include <WithoutReplacementResamplingAlgorithm.h>
#include <BestParticlesResamplingAlgorithm.h>
#include <FlatPowerProfile.h>
#include <ExponentialPowerProfile.h>
#include <MLSDmAlgorithm.h>
#include <PSPAlgorithm.h>
#include <TimeInvariantChannel.h>

class ISWCS10System : public ChannelOrderEstimationSystem
{
protected:
    int nSurvivors;
    bool adjustParticlesNumberFromSurvivors,adjustSurvivorsFromParticlesNumber;

	double _velocity;
	double _carrierFrequency;
	double _period;

	// vectors of channel estimators for the rows of the channel matrix (one channer per output algorithm)
	vector<ChannelMatrixEstimator *> kalmanChannelEstimators;
	
	// vectors of channel estimators for unknown channel order algorithms (single channel order)
	vector<ChannelMatrixEstimator *> kalmanWholeChannelEstimators;

    ResamplingAlgorithm *withoutReplacementResamplingAlgorithm,*bestParticlesResamplingAlgorithm;

    KalmanEstimator *_kalmanEstimatorForActualChannelOrder,*_kalmanEstimatorForMaximumChannelOrder;
	
	std::vector<uint> _subchannelOrders;

// 	// needed for solving ambiguity when computing MSE
// 	
// 	// each element of the vector (a vector itself) stores the signs of the best permutation corresponding
// 	// to one of the pieces into which the frame was dividing in order to compute the SER
// 	std::vector<std::vector<int> > _piecesBestPermutationSigns;
// 	
// 	// each element of the vector stores the index of the best permutation corresponding to one of the
// 	// pieces into which the frame was dividing in order to compute the SER
// 	std::vector<uint> _piecesBestPermuationIndexes;
// 	
// 	// a vector that stores the indexes within the frame (from 0 to "frameLength") where a sign change takes place
// 	// it ALWAYS includes the first index (0) and the last ("frameLength")
// 	std::vector<uint> _signChanges;

    virtual void addAlgorithms();
    virtual void buildSystemSpecificVariables();
    virtual void beforeEndingFrame();
public:
    ISWCS10System();

    ~ISWCS10System();
	
// 	virtual double computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,uint &iBestPermutation,vector<int> &bestPermutationSigns);
};

#endif
