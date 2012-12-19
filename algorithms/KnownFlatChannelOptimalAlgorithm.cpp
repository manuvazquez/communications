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
#include "KnownFlatChannelOptimalAlgorithm.h"

#include <assert.h>

KnownFlatChannelOptimalAlgorithm::KnownFlatChannelOptimalAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, const MIMOChannel& channel, uint preambleLength): KnownChannelAlgorithm(name, alphabet, L, Nr, N, iLastSymbolVectorToBeDetected, channel),_preambleLength(preambleLength),_detectedSymbols(_nInputs,iLastSymbolVectorToBeDetected-preambleLength)
{
    // a new alphabet extended with 0 (that meaning, no symbol is transmitted)
    vector<tSymbol> extendedAlphabetSymbols(_alphabet.length()+1);
    
    for(uint i=0;i<_alphabet.length();i++)
        extendedAlphabetSymbols[i] = _alphabet[i];
    extendedAlphabetSymbols[_alphabet.length()] = 0.0;
    
    _extendedAlphabet = new Alphabet(extendedAlphabetSymbols);
}

KnownFlatChannelOptimalAlgorithm::~KnownFlatChannelOptimalAlgorithm()
{
    delete _extendedAlphabet;
}

void KnownFlatChannelOptimalAlgorithm::run(MatrixXd observations, vector< double > noiseVariances)
{
    uint iAlphabet,i;
	uint childrenHeight;
    double UxS,newCost;
        
    for(uint iProcessedObservation=_preambleLength;iProcessedObservation<_iLastSymbolVectorToBeDetected;iProcessedObservation++)
    {
		// a new tree must be built for detecting the symbols transmitted at this time instant
		Tree tree;
		
		// the root node is initialized...
		KnownFlatChannelOptimalAlgorithmNode* currentNode = new KnownFlatChannelOptimalAlgorithmNode(_nInputs);
		
		// ...and added to the tree
		tree.setRoot(currentNode);
        
        // the corresponding channel matrix is kept in a variable (for the sake of clarity)
        MatrixXd H = _channel.getTransmissionMatrix(iProcessedObservation);
        
        // the Cholesky decomposition of HtH
        Eigen::LLT<MatrixXd> lltOfHTH(H.transpose()*H);
        
		VectorXd transformedObs = lltOfHTH.matrixL().solve(H.transpose()*observations.col(iProcessedObservation));
        
		MatrixXd U = lltOfHTH.matrixL().transpose();
        
        while(currentNode->getHeight()<_nInputs)
        {
			// all the children will have this height
			childrenHeight = currentNode->getHeight() + 1;
			
            for(iAlphabet=0;iAlphabet<getAlphabetAt(iProcessedObservation,childrenHeight)->length();iAlphabet++)
            {
				// the symbols vector from the current node is obtained...
				VectorXd newSymbolsVector = currentNode->getSymbolsVector();
				
				// ...and a new symbols is added to the end of it
				newSymbolsVector(_nInputs-childrenHeight) = getAlphabetAt(iProcessedObservation,childrenHeight)->operator[](iAlphabet);
                
                UxS = 0.0;
                for(i=0;i<childrenHeight;i++)
                    UxS += U(_nInputs-childrenHeight,_nInputs-1-i)*newSymbolsVector(_nInputs-1-i);
				
				// the cost of the corresponding child node is computed
				newCost = currentNode->getCost()+(transformedObs(_nInputs-childrenHeight)-UxS)*(transformedObs(_nInputs-childrenHeight)-UxS);
                                    
				// a child node is built with the updated symbols vector and cost
				KnownFlatChannelOptimalAlgorithmNode *child = new KnownFlatChannelOptimalAlgorithmNode(newCost,newSymbolsVector);
                
				// the child is added to the tree
				tree.addNode(child,currentNode);
            }
            
				// the best leaf node is chosen to continue
				currentNode = dynamic_cast<KnownFlatChannelOptimalAlgorithmNode *> (tree.bestLeaf());
        }
        
        VectorXd bestSymbolsVector = currentNode->getSymbolsVector();
		
		_detectedSymbols.col(iProcessedObservation) = bestSymbolsVector;
        
    } //for(uint iProcessedObservation=_preambleLength;iProcessedObservation<_iLastSymbolVectorToBeDetected;iProcessedObservation++)
}

MatrixXd KnownFlatChannelOptimalAlgorithm::getDetectedSymbolVectors()
{
    return _detectedSymbols;
}
