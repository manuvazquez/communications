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
#ifndef VITERBIALGORITHM_H
#define VITERBIALGORITHM_H

#include <KnownChannelAlgorithm.h>
#include <StillMemoryMIMOChannel.h>

/**
    @author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <Trellis.h>
#include <ViterbiPath.h>

enum tStage {exitStage,arrivalStage};

class ViterbiAlgorithm : public KnownChannelAlgorithm
{
protected:
    uint _d;
    Trellis *_trellis;
    ViterbiPath *_exitStage, *_arrivalStage;
    MatrixXd _preamble,*_detectedSymbolVectors;
	
	//! this is needed for algorithms that are children of Viterbi and require some kind of special initialization involving the first observation
	//! e.g. ViterbiAlgorithmWithAprioriProbabilities
	uint _iFirstInLoopProcessedObservation;

    virtual void deployState(int iState,const VectorXd &observations,const MatrixXd &channelMatrix,const double noiseVariance);
	void swapStages();
	
	/*!
	  it peforms the basic operations of any Viterbi algorithm (expansion of nodes, keep the best...). Specific operations for a particular Viterbi algorithm must be carried out in \ref run
	*/
	void process(MatrixXd observations,vector<double> noiseVariances,uint firstSymbolVectorDetectedAt);

public:
    ViterbiAlgorithm(string name, Alphabet alphabet,uint L,uint Nr,uint N, uint iLastSymbolVectorToBeDetected, const StillMemoryMIMOChannel& channel,const MatrixXd &preamble,uint smoothingLag);

    ~ViterbiAlgorithm();

    uint bestState()
    {
        uint bestState = 0;
        double bestCost = _exitStage[0].getCost();

        for(uint iState=1;iState<_trellis->nStates();iState++)
            if(_exitStage[iState].getCost() < bestCost)
            {
                bestState = iState;
                bestCost = _exitStage[iState].getCost();
            }
        return bestState;
    }
    
    void run(MatrixXd observations, vector< double > noiseVariances);    

    // detection will not start until the "firstSymbolVectorDetectedAt" observation
    virtual void run(MatrixXd observations,vector<double> noiseVariances,uint firstSymbolVectorDetectedAt);
    
    MatrixXd getDetectedSymbolVectors();
    void printStage(tStage exitOrArrival);
};

#endif
