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
#ifndef LINEARFILTERBASEDALGORITHM_H
#define LINEARFILTERBASEDALGORITHM_H

#include <KnownChannelOrderAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <LinearDetector.h>

class LinearFilterBasedAlgorithm : public KnownChannelOrderAlgorithm
{
protected:
	uint _d;
	LinearDetector *_linearDetector;
    MatrixXd _detectedSymbolVectors;
    std::vector<MatrixXd> _estimatedChannelMatrices;
	std::vector<double> _ARcoefficients;

    bool _substractContributionFromKnownSymbols;
	
	virtual void process(const MatrixXd &observations,vector<double> noiseVariances, MatrixXd trainingSequence);

public:
    LinearFilterBasedAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, LinearDetector *linearDetector, std::vector<double> ARcoefficients, bool substractContributionFromKnownSymbols=false);

    ~LinearFilterBasedAlgorithm();

    virtual void run(MatrixXd observations,vector<double> noiseVariances);
    virtual void run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence);

    virtual MatrixXd getDetectedSymbolVectors();
    
    virtual vector<MatrixXd> getEstimatedChannelMatrices();     

};

#endif
