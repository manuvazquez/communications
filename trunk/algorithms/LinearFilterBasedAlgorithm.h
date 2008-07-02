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
private:
	void Process(const tMatrix &observations,vector<double> noiseVariances, tMatrix trainingSequence);
protected:
	int _c,_d;
	LinearDetector *_linearDetector;
	tMatrix _detectedSymbolVectors;
	tMatrix *_estimatedChannelMatrices;
	double _ARcoefficient;

    bool _substractContributionFromKnownSymbols;

// 	virtual vector<tMatrix> EstimateChannelFromTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence);
public:
    LinearFilterBasedAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int backwardsSmoothingLag, int smoothingLag, LinearDetector *linearDetector, double ARcoefficient, bool substractContributionFromKnownSymbols=false);

    ~LinearFilterBasedAlgorithm();

    void Run(tMatrix observations,vector<double> noiseVariances);
    void Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);

    tMatrix GetDetectedSymbolVectors();
    vector<tMatrix> GetEstimatedChannelMatrices();

};

#endif
