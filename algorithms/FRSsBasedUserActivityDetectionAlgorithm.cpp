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


#include "FRSsBasedUserActivityDetectionAlgorithm.h"


#include <Eigen/QR>
#include <Octave.h>

// #define DEBUG

FRSsBasedUserActivityDetectionAlgorithm::FRSsBasedUserActivityDetectionAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, uint m, MatrixXd preamble, MatrixXd spreadingCodes, const std::vector<double> grid, const std::vector<UsersActivityDistribution> usersActivityPdfs, std::string channelTransitionProbabilitiesFileName): KnownChannelOrderAlgorithm(name, alphabet, L, Nr, N, iLastSymbolVectorToBeDetected,m,preamble),_spreadingCodes(spreadingCodes),_grid(grid),
_detectedSymbolVectors(N,iLastSymbolVectorToBeDetected),_estimatedChannelMatrices(iLastSymbolVectorToBeDetected),_estimatedChannelMatricesCells(iLastSymbolVectorToBeDetected,std::vector<uint>(N)),_usersActivityPdfs(usersActivityPdfs)
{
	// QR decomposition of the spreading codes matrix
	HouseholderQR<MatrixXd> qr(_spreadingCodes);
	_Qtrans = (qr.householderQ()).transpose();
	_R = qr.matrixQR().triangularView<Upper>();
	assert(_spreadingCodes.isApprox(_Qtrans.transpose()*_R));
	
#ifdef DEBUG
	cout << "spreading codes is" << endl << _spreadingCodes << endl;
	cout << "Q:" << endl << _Qtrans.transpose() << endl << "R:" << endl << _R << endl;
	cout << "Q' x Q:" << endl << _Qtrans*_Qtrans.transpose() << endl;
	cout << "the result of Q x R" << endl << _Qtrans.transpose()*_R << endl;
	getchar();
#endif
	
	_iFirstSymbolVectorToBeDetected = _preamble.cols();
	
	// the estimated channel transition probabilities are read from the received file name
	std::ifstream f;
	f.open(channelTransitionProbabilitiesFileName.c_str(),std::ifstream::in);
	if(f.fail())
		throw RuntimeException("FRSsBasedUserActivityDetectionAlgorithm::FRSsBasedUserActivityDetectionAlgorithm: error reading file \"" + channelTransitionProbabilitiesFileName + "\"");
	_estimatedChannelTransitionProbabilities = Octave::eigenFromOctaveFileStream(f);
	f.close();
	
	VectorXd sums = _estimatedChannelTransitionProbabilities.rowwise().sum();
	
	for(uint i=0;i<_estimatedChannelTransitionProbabilities.rows();i++)
		for(uint j=0;j<_estimatedChannelTransitionProbabilities.cols();j++)
			_estimatedChannelTransitionProbabilities(i,j) /= sums(i);
	
// 	cout << "estimatedChannelTransitionProbabilities" << endl << _estimatedChannelTransitionProbabilities << endl;
}

std::vector< MatrixXd> FRSsBasedUserActivityDetectionAlgorithm::getEstimatedChannelMatrices()
{
	return _estimatedChannelMatrices;
}

MatrixXd FRSsBasedUserActivityDetectionAlgorithm::getDetectedSymbolVectors()
{
	return _detectedSymbolVectors;
}

void FRSsBasedUserActivityDetectionAlgorithm::run(MatrixXd observations, std::vector< double> noiseVariances, MatrixXd trainingSequence)
{
	throw RuntimeException("FRSsBasedUserActivityDetectionAlgorithm::run: this is not implemented.");
}

void FRSsBasedUserActivityDetectionAlgorithm::run(MatrixXd observations, std::vector< double> noiseVariances)
{
	uint i,iAlphabet,iCell,iCurrentUser,childrenHeight;
	double RxD,No,newCost;
	VectorXd newSymbolsVector;
	MatrixXd newChannelMatrix;
	std::vector<uint> newChannelMatrixCells;
	
	Alphabet extendedAlphabet = _alphabet.buildNewAlphabetByAddingSymbol(0.0);
	
	for(uint iObservationToBeProcessed=_iFirstSymbolVectorToBeDetected;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
	{
#ifdef DEBUG
		cout << "iObservationToBeProcessed = " << iObservationToBeProcessed << endl;
#endif
		No = noiseVariances[iObservationToBeProcessed]*2;
		
		// the "transformed" observation is first computed using the QR decomposition
		VectorXd transformedObs = _Qtrans*observations.col(iObservationToBeProcessed);
		
		Tree tree;
		
		// the root node is initialized...
		FRSsBasedUserActivityDetectionAlgorithmNode *currentNode = new FRSsBasedUserActivityDetectionAlgorithmNode(_nInputs);
		
		// ...and added to the tree
		tree.setRoot(currentNode);

		while(currentNode->getHeight()<_nInputs)
        {
			// the height of the children nodes is one unit bigger than that of the parent
			childrenHeight = currentNode->getHeight() + 1;
			
			// for the sake of clarity...
			iCurrentUser = _nInputs-childrenHeight;
			
            for(iAlphabet=0;iAlphabet<extendedAlphabet.length();iAlphabet++)
            {
				// ...also for the sake of clarity
				double currentSymbol = extendedAlphabet[iAlphabet];
				
				// the symbols vector for the children of the current node is obtained from the latter...
				newSymbolsVector = currentNode->getSymbolsVector();

				// ...and updated with the symbol being tested
				newSymbolsVector(iCurrentUser) = currentSymbol;
				
				for(iCell=0;iCell<_grid.size();iCell++)
				{	
					// the channel matrix coefficients (and their associated cells) of the children are obtained from those of the parent node (the current node)...
					newChannelMatrix = currentNode->getChannelMatrix();
					newChannelMatrixCells = currentNode->getChannelMatrixCells();
					
					
					// ...and updated
					newChannelMatrix(iCurrentUser) = _grid[iCell];
					newChannelMatrixCells[iCurrentUser] = iCell;
					
					RxD = 0.0;
					for(i=0;i<childrenHeight;i++)
						RxD += _R(iCurrentUser,_nInputs-1-i)*newSymbolsVector(_nInputs-1-i)*newChannelMatrix(_nInputs-1-i);
				
					// the new cost is also initialized to that of the parent
					newCost = currentNode->getCost();
						
					newCost += (transformedObs(iCurrentUser)-RxD)*(transformedObs(iCurrentUser)-RxD);
					
					if(iObservationToBeProcessed==_iFirstSymbolVectorToBeDetected)
						newCost += -No*log(_usersActivityPdfs[iCurrentUser].probApriori(Util::isUserActive(currentSymbol)));
					else
						newCost += -No*log(_usersActivityPdfs[iCurrentUser].probXgivenY(Util::isUserActive(currentSymbol),Util::isUserActive(_detectedSymbolVectors(iCurrentUser,iObservationToBeProcessed-1))));
					
					newCost += No*double(Util::isUserActive(currentSymbol))*log(_alphabet.length());
					
					// if the user is active at the corresponding time instant...
					if(Util::isUserActive(currentSymbol))
					{
						// ...and this is the first time instant or the user was not active before...
						if(iObservationToBeProcessed==_iFirstSymbolVectorToBeDetected || !Util::isUserActive(_detectedSymbolVectors(iCurrentUser,iObservationToBeProcessed-1)))
							// ...then the probability of its corresponding channel coefficient is given by the "a priori" probability of a channel coefficient
							newCost += -No*log(channelCoeffAprioriProb(newChannelMatrix(iCurrentUser)));
						// ...otherwise it's a surviving user
						else
						{
							// ...and the probability of the channel coefficient at the current time instant is conditional on the value it had before
							newCost += -No*log(channelCoeffConditionalProb(newChannelMatrixCells[iCurrentUser],_estimatedChannelMatricesCells[iObservationToBeProcessed-1][iCurrentUser]));
						}
					}
					
#ifdef DEBUG
					cout << "----------------- " << "iCurrentUser = " << iCurrentUser << " (child height = " << childrenHeight << ") [symbol = " << child.symbolsVector(iCurrentUser) << ", channel coeff. = " << child.channelMatrix(iCurrentUser) << "]" << endl;
					cout << "transformedObs(iCurrentUser) = " << transformedObs(iCurrentUser) << ", RxD = " << RxD << ", child.cost = " << child.cost << endl;
#endif
					
					FRSsBasedUserActivityDetectionAlgorithmNode *child = new FRSsBasedUserActivityDetectionAlgorithmNode(newCost,newSymbolsVector,newChannelMatrix,newChannelMatrixCells);
					
					// the child is added to the tree
					tree.addNode(child,currentNode);
				}
            } // for(iAlphabet=0;iAlphabet<extendedAlphabet.length();iAlphabet++)
			
			// the best leaf node is chosen to continue
			currentNode = dynamic_cast<FRSsBasedUserActivityDetectionAlgorithmNode *> (tree.bestLeaf());
			
#ifdef DEBUG
			cout << "======= iCurrentUser = " << iCurrentUser << " ==========" << endl;
			cout << "best node has symbols:" << endl << nodes[iCurrentNode].symbolsVector << endl << " and channel coeffs: " << endl << nodes[iCurrentNode].channelMatrix << endl;
			getchar();
#endif
        } // while(currentNode->getHeight()<_nInputs)

		_detectedSymbolVectors.col(iObservationToBeProcessed) = currentNode->getSymbolsVector();
		_estimatedChannelMatrices[iObservationToBeProcessed] = currentNode->getChannelMatrix();
		_estimatedChannelMatricesCells[iObservationToBeProcessed] = currentNode->getChannelMatrixCells();

#ifdef DEBUG
		cout << "transformed observations:" << endl << transformedObs << endl;
		cout << "detected symbols: " << endl << _detectedSymbolVectors.col(iObservationToBeProcessed) << endl;
		cout << "detected channel matrix: " << endl << _estimatedChannelMatrices[iObservationToBeProcessed] << endl;
		cout << "-------------------------- iObservationToBeProcessed = " << iObservationToBeProcessed << " ----------------------------- " << endl;
		getchar();
		
#endif

	} // for(uint iObservationToBeProcessed=_iFirstSymbolVectorToBeDetected;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

double FRSsBasedUserActivityDetectionAlgorithm::channelCoeffAprioriProb(double channelCoeff)
{
	return 1.0/double(_grid.size());
}

double FRSsBasedUserActivityDetectionAlgorithm::channelCoeffConditionalProb(uint currentChannelCoeffCell, uint previousChannelCoeffCell)
{
	return _estimatedChannelTransitionProbabilities(previousChannelCoeffCell,currentChannelCoeffCell);
}