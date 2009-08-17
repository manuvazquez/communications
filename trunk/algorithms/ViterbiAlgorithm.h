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
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
// #include <lapackpp/blas3pp.h>

enum tStage {exitStage,arrivalStage};

class ViterbiAlgorithm : public KnownChannelAlgorithm
{
private:
	// decimal inputs will be converted to a symbol vector and stored in here
	vector<tSymbol> _inputVector;

	// states in decimal format will be converted to a symbol vector and stored in here
	vector<tSymbol> _stateVector;
protected:
    int _d;
	Trellis _trellis;
    ViterbiPath *_exitStage, *_arrivalStage;
    MatrixXd _preamble,*_detectedSymbolVectors;

    void DeployState(int iState,const VectorXd &observations,const MatrixXd &channelMatrix);
public:
    ViterbiAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected, const StillMemoryMIMOChannel& channel,const tMatrix &preamble,int smoothingLag);

    ~ViterbiAlgorithm();

    int BestState()
    {
        int bestState = 0;
        double bestCost = _exitStage[0].getCost();

        for(int iState=1;iState<_trellis.Nstates();iState++)
            if(_exitStage[iState].getCost() < bestCost)
            {
                bestState = iState;
                bestCost = _exitStage[iState].getCost();
            }
        return bestState;
    }
    
    void run(tMatrix observations, vector< double > noiseVariances)
    {
        run(Util::lapack2eigen(observations),noiseVariances);
    }
    void run(MatrixXd observations, vector< double > noiseVariances);    

    // detection will not start until the "firstSymbolVectorDetectedAt" observation
    void run(MatrixXd observations,vector<double> noiseVariances,int firstSymbolVectorDetectedAt);
    void run(tMatrix observations, vector< double > noiseVariances,int firstSymbolVectorDetectedAt)
    {
        run(Util::lapack2eigen(observations),noiseVariances,firstSymbolVectorDetectedAt);
    }    
    
    tMatrix getDetectedSymbolVectors()
    {
        return Util::eigen2lapack(getDetectedSymbolVectors_eigen());
    }
    MatrixXd getDetectedSymbolVectors_eigen();
    void PrintStage(tStage exitOrArrival);
};

#endif
