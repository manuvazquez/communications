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
#ifndef CMEBASEDALGORITHM_H
#define CMEBASEDALGORITHM_H

#include <UnknownChannelOrderAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <RLSEstimator.h>

class CMEBasedAlgorithm : public UnknownChannelOrderAlgorithm
{
protected:
    tMatrix _symbolVectors;
public:
    CMEBasedAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, const tMatrix &symbolVectors);

	virtual void Run(tMatrix observations,vector<double> noiseVariances);
	virtual void Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);

    virtual tMatrix GetDetectedSymbolVectors();
    virtual vector<tMatrix> GetEstimatedChannelMatrices();

//     bool PerformsChannelOrderAPPEstimation() const { return false;}

    ~CMEBasedAlgorithm();

};

#endif
