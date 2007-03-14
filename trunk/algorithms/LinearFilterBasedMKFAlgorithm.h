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
#ifndef LINEARFILTERBASEDMKFALGORITHM_H
#define LINEARFILTERBASEDMKFALGORITHM_H

#include <LinearFilterBasedSMCAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <KalmanEstimator.h>

class LinearFilterBasedMKFAlgorithm : public LinearFilterBasedSMCAlgorithm
{
public:
    LinearFilterBasedMKFAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, KalmanEstimator* channelEstimator, LinearDetector* linearDetector, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances, double ARcoefficient, double samplingVariance, double ARprocessVariance);

protected:
    virtual void FillFirstEstimatedChannelMatrix(int iParticle, tMatrix& firstEstimatedChannelMatrix)
    {
    	firstEstimatedChannelMatrix = (dynamic_cast<KalmanEstimator *> (_particleFilter->GetParticle(iParticle)->GetChannelMatrixEstimator(_estimatorIndex)))->SampleFromPredictive();
    }

};

#endif
