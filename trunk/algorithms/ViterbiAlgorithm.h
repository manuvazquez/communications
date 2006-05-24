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

/**
	@author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
// #include <lapackpp/blas3pp.h>

enum tStage {exitStage,arrivalStage};

class ViterbiAlgorithm : public KnownChannelAlgorithm
{
protected:
	int _nStates,_nPossibleInputs;
	int **_stateTransitionMatrix;

	typedef struct {
		double cost;
		tMatrix *sequence;
	} tState;

	tState *_exitStage, *_arrivalStage;
	tMatrix _preamble;
	tRange rAllSymbolRows,rmMinus1FirstColumns;

	void BuildStateTransitionMatrix();
	void DeployState(int iState,const tVector &observations,const tMatrix &channelMatrix);
public:
    ViterbiAlgorithm(string name, Alphabet alphabet, const MIMOChannel& channel,const tMatrix &preamble);

    ~ViterbiAlgorithm();

	int BestState()
	{
		int bestState = 0;
		double bestCost = _exitStage[0].cost;
	
		for(int iState=1;iState<_nStates;iState++)
			if(_exitStage[iState].cost < bestCost)
			{
				bestState = iState;
				bestCost = _exitStage[iState].cost;
			}
		return bestState;
	}
	void Run(const tMatrix &observations,vector<double> noiseVariances);
	void Run(const tMatrix &observations,vector<double> noiseVariances,int detectionLag);
	void PrintStage(tStage exitOrArrival);
	double SER(tMatrix symbols);
};

#endif
