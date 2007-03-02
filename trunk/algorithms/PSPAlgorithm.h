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
#ifndef PSPALGORITHM_H
#define PSPALGORITHM_H

#include <KnownChannelOrderAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <Trellis.h>
#include <PSPPath.h>
#include <PSPPathCandidate.h>

class PSPAlgorithm : public KnownChannelOrderAlgorithm
{
private:
	tRange _rAllSymbolRows;

	// decimal inputs will be converted to a symbol vector and stored in here
	vector<tSymbol> _inputVector;

	// states in decimal format will be converted to a symbol vector and stored in here
	vector<tSymbol> _stateVector;

	void BestPairStateSurvivor(int &bestState,int &bestSurvivor);
protected:
    int _nSurvivors,_d,_startDetectionTime;
	Trellis _trellis;
    PSPPath **_exitStage, **_arrivalStage;
    tMatrix *_detectedSymbolVectors;
    std::vector<tMatrix> _estimatedChannelMatrices;
	int _firstSymbolVectorDetectedAt;
	double _ARcoefficient;
	PSPPathCandidate **_bestArrivingPaths;

	void ProcessOneObservation(const tVector &observations,double noiseVariance);
	void Process(const tMatrix &observations,vector<double> noiseVariances);
	void DeployState(int iState,const tVector &observations, double noiseVariance);
public:
    PSPAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int firstSymbolVectorDetectedAt, double ARcoefficient);

    ~PSPAlgorithm();

	void Run(tMatrix observations,vector<double> noiseVariances);
	void Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);

	tMatrix GetDetectedSymbolVectors();
	std::vector<tMatrix> GetEstimatedChannelMatrices();

};

#endif
