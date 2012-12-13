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

// #define DEBUG

FRSsBasedUserActivityDetectionAlgorithm::FRSsBasedUserActivityDetectionAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, uint m, MatrixXd preamble, MatrixXd spreadingCodes, double firstCell, double lastCell, uint nCells, const std::vector<UsersActivityDistribution> usersActivityPdfs): KnownChannelOrderAlgorithm(name, alphabet, L, Nr, N, iLastSymbolVectorToBeDetected,m,preamble),_spreadingCodes(spreadingCodes),
_detectedSymbolVectors(N,iLastSymbolVectorToBeDetected),_estimatedChannelMatrices(iLastSymbolVectorToBeDetected),_usersActivityPdfs(usersActivityPdfs)
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
#endif
	
	_iFirstSymbolVectorToBeDetected = _preamble.cols();
	
	_grid = std::vector<double>(nCells);
	_gridStep = (lastCell-firstCell)/(nCells-1);
	
	_grid[0] = firstCell;
	for(uint i=1;i<(nCells-1);i++)
		_grid[i] = _grid[i-1] + _gridStep;
	_grid[nCells-1] = lastCell;
	
#ifdef DEBUG
	cout << "grid step = " << _gridStep << endl;
	cout << "grid:" << endl << _grid << endl;
#endif
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
	uint i,iAlphabet,iCell,iCurrentNode,childrenHeight;
	double RxD,No;
	
	Alphabet extendedAlphabet = _alphabet.buildNewAlphabetByAddingSymbol(0.0);
	
	for(uint iObservationToBeProcessed=_iFirstSymbolVectorToBeDetected;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
	{
#ifdef DEBUG
		cout << "iObservationToBeProcessed = " << iObservationToBeProcessed << endl;
#endif
		No = noiseVariances[iObservationToBeProcessed]*2;
		
		// the "transformed" observation is first computed using the QR decomposition
		VectorXd transformedObs = _Qtrans*observations.col(iObservationToBeProcessed);
		
		std::vector<tTreeNode> nodes;
		
        // root node is initialized
        tTreeNode rootNode;
        rootNode.cost = 0.0;
        rootNode.height = 0;
        rootNode.id = 0;
		rootNode.symbolsVector = VectorXd::Zero(_nInputs);
		rootNode.channelMatrix = MatrixXd::Zero(_Nr,_nInputs);
		
        // root node is added to the list of nodes
        nodes.push_back(rootNode);

        // we start by the root node
        iCurrentNode = 0;    
        
        while(nodes[iCurrentNode].height<_nInputs)
        {
			// the height of the children nodes is one unit bigger than that of the parent
            childrenHeight = nodes[iCurrentNode].height+1;
			
			// for the sake of clarity...
			uint iCurrentUser = _nInputs-childrenHeight;
			
            for(iAlphabet=0;iAlphabet<extendedAlphabet.length();iAlphabet++)
            {
				// ...also for the sake of clarity
				double currentSymbol = extendedAlphabet[iAlphabet];
				
				for(iCell=0;iCell<_grid.size();iCell++)
				{
#ifdef DEBUG
// 					cout << "iCurrentUser = " << iCurrentUser << endl;
// 					cout << "extendedAlphabet[iAlphabet] = " << extendedAlphabet[iAlphabet] << " _grid[iCell] = " << _grid[iCell] << endl;
#endif
					// the parent node is replicated
					tTreeNode child = nodes[iCurrentNode];
					
					child.height = childrenHeight;
					child.symbolsVector(iCurrentUser) = currentSymbol;
					child.channelMatrix(iCurrentUser) = _grid[iCell];
					
					RxD = 0.0;
				
					for(i=0;i<child.height;i++)
						RxD += _R(iCurrentUser,_nInputs-1-i)*child.symbolsVector(_nInputs-1-i)*child.channelMatrix(_nInputs-1-i);
										
					child.cost += (transformedObs(iCurrentUser)-RxD)*(transformedObs(iCurrentUser)-RxD);
					
					if(iObservationToBeProcessed==_iFirstSymbolVectorToBeDetected)
					{
// 						child.cost += -No*double(Util::isUserActive(currentSymbol))*log(_usersActivityPdfs[iCurrentUser].probApriori(true));
						child.cost += -No*log(_usersActivityPdfs[iCurrentUser].probApriori(Util::isUserActive(currentSymbol)));
// 						child.cost += No*double(Util::isUserActive(currentSymbol))*log(_alphabet.length());
					}else
						child.cost += -No*log(_usersActivityPdfs[iCurrentUser].probXgivenY(Util::isUserActive(currentSymbol),Util::isUserActive(_detectedSymbolVectors(iCurrentUser,iObservationToBeProcessed-1))));
					
					child.cost += No*double(Util::isUserActive(currentSymbol))*log(_alphabet.length());
					
					// if the user is active at the corresponding time instant...
					if(Util::isUserActive(currentSymbol))
					{
						// ...and this is the first time instant or the user was not active before...
						if(iObservationToBeProcessed==_iFirstSymbolVectorToBeDetected || !Util::isUserActive(_detectedSymbolVectors(iCurrentUser,iObservationToBeProcessed-1)))
							// ...then the probability of its corresponding channel coefficient is given by the "a priori" probability of a channel coefficient
							child.cost += -No*log(channelCoeffAprioriProb(child.channelMatrix(iCurrentUser)));
						// ...otherwise it's a surviving user
						else
						{
							// ...and the probability of the channel coefficient at the current time instant is conditional on the value it had before
							child.cost += -No*log(channelCoeffConditionalProb(child.channelMatrix(iCurrentUser),_estimatedChannelMatrices[iObservationToBeProcessed-1](iCurrentUser)));
						}
					}
							
					
					// the children vector of this node is cleared (the one inherited from the father because of the copy may not be empty if this is not the the first child
					child.children.clear();
					
					// child will be in this position of the vector
					child.id = nodes.size();
					
					// parent node is updated
					nodes[iCurrentNode].children.push_back(child.id);
					
					// node is added to the list
					nodes.push_back(child);
				}
            }
            
            // best node is chosen
            iCurrentNode = iBestLeaf(nodes);
        }
        
		_detectedSymbolVectors.col(iObservationToBeProcessed) = nodes[iCurrentNode].symbolsVector;
		_estimatedChannelMatrices[iObservationToBeProcessed] = nodes[iCurrentNode].channelMatrix;
        
        // for the next iteration
        nodes.clear();
		
#ifdef DEBUG
		cout << "-----------------------------" << endl;
#endif

	} // for(uint iObservationToBeProcessed=_iFirstSymbolVectorToBeDetected;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

uint FRSsBasedUserActivityDetectionAlgorithm::iBestLeaf(const vector<tTreeNode> &nodes)
{
	assert(nodes.size()>0);
    uint iBest;
    double bestCost = 0.0;

	uint i=0;
	
	// the first leaf node is searched for
	while(i<nodes.size() && nodes[i].children.size()!=0)
		i++;
	
	iBest = i;
	bestCost = nodes[i].cost;
	
    for(i=i+1;i<nodes.size();i++)
    {
        // if it isn't a leaf node...
        if(nodes[i].children.size()!=0)
            continue;

        if(nodes[i].cost < bestCost)
        {
            iBest = i;
            bestCost = nodes[i].cost;
        }
    }
        
    return iBest;
}

