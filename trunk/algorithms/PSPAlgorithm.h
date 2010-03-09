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
	// decimal inputs will be converted to a symbol vector and stored in here
	vector<tSymbol> _inputVector;

	// states in decimal format will be converted to a symbol vector and stored in here
	vector<tSymbol> _stateVector;

	/*!
	  It goes through all the survivors in all the states of \ref _exitStage and finds the one with the smaller cost
	  \param bestState returned
	  \param bestSurvivor returned
	*/
	void bestPairStateSurvivor(int &bestState,int &bestSurvivor);
	
	/*!
	  It seeks the survivor with the worst cost...or simply an empty slot in the list of \ref _bestArrivingPaths
	  \param iState the state of \ref _bestArrivingPaths that is to be searched in
	  \return the index of the survivor with the worst cost
	*/
    int disposableSurvivor(int iState);
	
	/*!
	  It contains code shared by the two \ref run methods
	*/
	void initializeTrellis();
protected:
    int _nSurvivors,_d,_startDetectionTime;
	Trellis *_trellis;
    PSPPath **_exitStage, **_arrivalStage;
    MatrixXd *_detectedSymbolVectors;
    std::vector<MatrixXd> _estimatedChannelMatrices;
	int _firstSymbolVectorDetectedAt;
	double _ARcoefficient;
	PSPPathCandidate **_bestArrivingPaths;

    void processOneObservation(const VectorXd &observations,double noiseVariance);
    void process(const MatrixXd &observations,vector<double> noiseVariances);
    virtual void deployState(int iState,const VectorXd &observations, double noiseVariance);
public:
    PSPAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, int smoothingLag, int firstSymbolVectorDetectedAt, double ARcoefficient, int nSurvivors);

    ~PSPAlgorithm();

    virtual void run(MatrixXd observations, vector< double > noiseVariances);
    virtual void run(MatrixXd observations, vector< double > noiseVariances, MatrixXd trainingSequence);   

    MatrixXd getDetectedSymbolVectors();
    std::vector<MatrixXd> getEstimatedChannelMatrices();
    void printState(int iState);

};

#endif
