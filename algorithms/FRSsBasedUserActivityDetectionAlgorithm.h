/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2012  Manu <manuavazquez@gmail.com>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#ifndef FRSSBASEDUSERACTIVITYDETECTIONALGORITHM_H
#define FRSSBASEDUSERACTIVITYDETECTIONALGORITHM_H

#include <KnownChannelOrderAlgorithm.h>

#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>

#include <UsersActivityDistribution.h>


class FRSsBasedUserActivityDetectionAlgorithm : public KnownChannelOrderAlgorithm
{
protected:
	const MatrixXd _spreadingCodes;
	MatrixXd _Qtrans,_R;
	
	const std::vector<double> _grid;
	
	uint _iFirstSymbolVectorToBeDetected;
	
    MatrixXd _detectedSymbolVectors;
	std::vector<MatrixXd> _estimatedChannelMatrices;
	std::vector<std::vector<uint> > _estimatedChannelMatricesCells;
	
	const std::vector<UsersActivityDistribution> _usersActivityPdfs; /// objects describing the pdf of the users activity
	
	MatrixXd _estimatedChannelTransitionProbabilities;
	
    typedef struct{
        double cost;
        std::vector<uint> children;
        uint height,id;
        VectorXd symbolsVector;
		MatrixXd channelMatrix;
		std::vector<uint> channelMatrixCells;
    } tTreeNode;
	
	uint iBestLeaf(const std::vector< FRSsBasedUserActivityDetectionAlgorithm::tTreeNode >& nodes);
	
// 	double channelCoeffAprioriProb(double channelCoeff) { return 1.0; }
// 	double channelCoeffConditionalProb(double currentChannelCoeff, double previousChannelCoeff) { return 1.0; }
	double channelCoeffAprioriProb(uint channelCoeff);
	double channelCoeffConditionalProb(uint currentChannelCoeffCell, uint previousChannelCoeffCell);
public:
    FRSsBasedUserActivityDetectionAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, uint m, MatrixXd preamble, MatrixXd spreadingCodes, const std::vector<double> grid, const std::vector<UsersActivityDistribution> usersActivityPdfs, MatrixXd estimatedChannelTransitionProbabilities);
    virtual std::vector< MatrixXd> getEstimatedChannelMatrices();
    virtual MatrixXd getDetectedSymbolVectors();
    virtual void run(MatrixXd observations, std::vector< double> noiseVariances, MatrixXd trainingSequence);
    virtual void run(MatrixXd observations, std::vector< double> noiseVariances);
};

#endif // FRSSBASEDUSERACTIVITYDETECTIONALGORITHM_H
