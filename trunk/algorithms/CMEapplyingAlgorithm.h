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
#ifndef CMEAPPLYINGALGORITHM_H
#define CMEAPPLYINGALGORITHM_H

#include <UnknownChannelOrderAlgorithm.h>

/**
    @author Manu <manu@rustneversleeps>
*/

#include <math.h>

#include <RLSEstimator.h>
#include <Algorithm.h>
#include <LinearFilterBasedAlgorithm.h>

class CMEapplyingAlgorithm : public UnknownChannelOrderAlgorithm
{
protected:
    std::vector<Algorithm *> algorithms;

    virtual tMatrix detectedSymbolsForChannelOrder(uint iChannelOrder,const tMatrix &observations,const vector<double> &noiseVariances,const tMatrix& trainingSequence) = 0;
    virtual std::vector<tMatrix> estimatedChannelMatricesForChannelOrder(uint iChannelOrder,const tMatrix &observations,const vector<double> &noiseVariances,const tMatrix& trainingSequence) = 0;
public:
    CMEapplyingAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected,vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble);

    virtual void run(tMatrix observations,vector<double> noiseVariances);
    virtual void run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);

    virtual tMatrix getDetectedSymbolVectors();
    virtual vector<tMatrix> getEstimatedChannelMatrices();
};

#endif
