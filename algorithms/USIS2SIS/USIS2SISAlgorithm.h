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
#ifndef USIS2SISALGORITHM_H
#define USIS2SISALGORITHM_H

#include <USIS.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <LinearFilterBasedSMCAlgorithm.h>
#include <TransitionCriterion.h>

class USIS2SISAlgorithm : public USIS
{
protected:
    TransitionCriterion *_transitionCriterion;
public:
    USIS2SISAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, vector< LinearDetector * > linearDetectors, MatrixXd preamble, uint iFirstObservation, uint smoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm, ChannelOrderEstimator* channelOrderEstimator, double ARcoefficient, double samplingVariance, double ARprocessVariance, TransitionCriterion *transitionCriterion);

protected:
    virtual void beforeResamplingProcess(uint iProcessedObservation, const MatrixXd& observations, const vector<double> &noiseVariances);

};

#endif
