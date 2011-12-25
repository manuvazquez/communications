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
#include <KalmanEstimator.h>

class PSPAlgorithm : public KnownChannelOrderAlgorithm
{
private:
	// decimal inputs will be converted to a symbol vector and stored in here
	vector<tSymbol> _inputVector;

	// states in decimal format will be converted to a symbol vector and stored in here
	vector<tSymbol> _stateVector;

	/*!
	  It goes through all the survivors in all the states of \ref _exitStage and finds the one with the smaller cost
	  \param bestState returned
	  \param bestSurvivor returned
	*/
	void bestPairStateSurvivor(uint &bestState,uint &bestSurvivor);
	
	/*!
	  It contains code shared by the two \ref run methods
	*/
	void initializeTrellis();
protected:
    uint _nSurvivors,_d,_startDetectionTime;
	Trellis *_trellis;
    PSPPath **_exitStage, **_arrivalStage;
    MatrixXd *_detectedSymbolVectors;
    std::vector<MatrixXd> _estimatedChannelMatrices;
	uint _firstSymbolVectorDetectedAt;
	PSPPathCandidate **_bestArrivingPaths;

	//! this variable will always be equal to \ref _startDetectionTime in this algorithm but not in \ref PSPAlgorithmWithAprioriProbabilities
	uint _iFirstInLoopProcessedObservation;

	/*!
	  It seeks the survivor with the worst cost...or simply an empty slot in the list of \ref _bestArrivingPaths
	  \param iState the state of \ref _bestArrivingPaths that is to be searched in
	  \return the index of the survivor with the worst cost
	*/
    uint disposableSurvivor(uint iState);

    void processOneObservation(const VectorXd &observations,double noiseVariance);
    void process(const MatrixXd &observations,vector<double> noiseVariances);
    virtual void deployState(uint iState,const VectorXd &observations, double noiseVariance);
public:
    PSPAlgorithm(string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, uint firstSymbolVectorDetectedAt, uint nSurvivors);

    ~PSPAlgorithm();

    virtual void run(MatrixXd observations, vector< double > noiseVariances);
    virtual void run(MatrixXd observations, vector< double > noiseVariances, MatrixXd trainingSequence);   

    MatrixXd getDetectedSymbolVectors();
    std::vector<MatrixXd> getEstimatedChannelMatrices();
    void printState(int iState);

};

#endif
