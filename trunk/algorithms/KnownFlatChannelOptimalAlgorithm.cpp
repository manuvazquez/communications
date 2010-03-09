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

// #define DEBUG3

KnownFlatChannelOptimalAlgorithm::KnownFlatChannelOptimalAlgorithm(string name, Alphabet alphabet, int L, int Nr, int N, int iLastSymbolVectorToBeDetected, const MIMOChannel& channel, int preambleLength): KnownChannelAlgorithm(name, alphabet, L, Nr, N, iLastSymbolVectorToBeDetected, channel),_preambleLength(preambleLength),_detectedSymbols(_nInputs,iLastSymbolVectorToBeDetected-preambleLength)
{
    // a new alphabet extended with 0 (that meaning, no symbol is transmitted)
    vector<tSymbol> extendedAlphabetSymbols(_alphabet.length()+1);
    
    for(int i=0;i<_alphabet.length();i++)
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
    int iAlphabet,iCurrentNode,i,childrenHeight;
    double UxS;
        
    vector<tTreeNode> nodes;
    
    for(int iProcessedObservation=_preambleLength;iProcessedObservation<_iLastSymbolVectorToBeDetected;iProcessedObservation++)
    {
        // root node is initialized
        tTreeNode rootNode;
        rootNode.cost = 0.0;
        
        rootNode.height = 0;
        rootNode.id = 0;
        rootNode.symbolsVector = VectorXd(_nInputs);
        
        // root node is added to the list of nodes
        nodes.push_back(rootNode);
        
        // the corresponding channel matrix is kept in a variable (for the sake of clarity)
        MatrixXd H = _channel.getTransmissionMatrix(iProcessedObservation);
        
        // the Cholesky decomposition of HtH
        Eigen::LLT<MatrixXd> lltOfHTH(H.transpose()*H);
        
        VectorXd transformedObs = lltOfHTH.matrixL().inverse()*H.transpose()*observations.col(iProcessedObservation);
        
        MatrixXd U = lltOfHTH.matrixL().transpose();
        
        // we start by the root node
        iCurrentNode = 0;    
        
        while(nodes[iCurrentNode].height<_nInputs)
        {
            childrenHeight = nodes[iCurrentNode].height+1;
            for(iAlphabet=0;iAlphabet<getAlphabetAt(iProcessedObservation,childrenHeight)->length();iAlphabet++)
            {
                // the parent node is replicated
                tTreeNode child = nodes[iCurrentNode];
                
                child.height = childrenHeight;
                child.symbolsVector(_nInputs-child.height) = (*getAlphabetAt(iProcessedObservation,childrenHeight))[iAlphabet];
                
                UxS = 0.0;
            
                for(i=0;i<child.height;i++)
                    UxS += U(_nInputs-child.height,_nInputs-1-i)*child.symbolsVector(_nInputs-1-i);
                                    
                child.cost = child.cost + (transformedObs(_nInputs-child.height)-UxS)*(transformedObs(_nInputs-child.height)-UxS);
                
                // the children vector of this node is cleared (the one inherited from the father because of the copy may not be empty if this is not the
                // the first child
                child.children.clear();
                
                // child will be in this position of the vector
                child.id = nodes.size();
                
                // parent node is updated
                nodes[iCurrentNode].children.push_back(child.id);
                
                // node is added to the list
                nodes.push_back(child);
            }
            
            // best node is chosen
            iCurrentNode = iBestLeaf(nodes);
        }
        
        for(i=0;i<_nInputs;i++)
            _detectedSymbols(i,iProcessedObservation-_preambleLength) = nodes[iCurrentNode].symbolsVector(i);
        
        // for the next iteration
        nodes.clear();
        
    } //for(int iProcessedObservation=_preambleLength;iProcessedObservation<_iLastSymbolVectorToBeDetected;iProcessedObservation++)
}

int KnownFlatChannelOptimalAlgorithm::iBestLeaf(const vector<tTreeNode> &nodes)
{
    int iBest = -1;
    double bestCost = 0.0;

    for(uint i=0;i<nodes.size();i++)
    {
        // if it isn't a leaf node
        if(nodes[i].children.size()!=0)
            continue;

        if(iBest==-1 || nodes[i].cost < bestCost)
        {
            iBest = i;
            bestCost = nodes[i].cost;
        }
    }
    
    return iBest;
}

MatrixXd KnownFlatChannelOptimalAlgorithm::getDetectedSymbolVectors()
{
    return _detectedSymbols;
}
