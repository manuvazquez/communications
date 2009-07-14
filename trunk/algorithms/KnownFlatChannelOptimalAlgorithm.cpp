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
}


void KnownFlatChannelOptimalAlgorithm::Run(tMatrix observations, vector< double > noiseVariances)
{
#ifdef DEBUG2
    cout << "en KnownFlatChannelOptimalAlgorithm::Run" << endl;
    cout << "observations are" << endl << observations;
#endif
    int iAlphabet,iCurrentNode,i;
    tMatrix HtH(_nInputs,_nInputs),invL_Ht(_nInputs,_nOutputs),U(_nInputs,_nInputs);
    tVector transformedObs(_nInputs);
    tLongIntVector piv(_nOutputs);
    double UxS;
        
    vector<tTreeNode> nodes;
    
    for(int iProcessedObservation=_preambleLength;iProcessedObservation<_iLastSymbolVectorToBeDetected;iProcessedObservation++)
    {
    
#ifdef DEBUG3
        cout << "--------- iProcessedObservation = " << iProcessedObservation << " ------------" << endl;
#endif    
        // root node is initialized
        tTreeNode rootNode;
        rootNode.cost = 0.0;
        
        rootNode.height = 0;
        rootNode.id = 0;
        rootNode.symbolsVector = tVector(_nInputs);
        
        // root node is added to the list of nodes
        nodes.push_back(rootNode);
        
        // the corresponding channel matrix is kept in a variable (for the sake of clarity)
        tMatrix H;
        H = _channel.getTransmissionMatrix(iProcessedObservation);
        
        // HtH = H'*H
        Blas_Mat_Trans_Mat_Mult(H,H,HtH);
        
        // the Cholesky decomposition of HtH
        tMatrix L = Util::cholesky(HtH);
        
        Util::transpose(L,U);
        
        // invL = inverse(L)
        tMatrix invL = L;
        LUFactorizeIP(invL,piv);
        LaLUInverseIP(invL,piv);
        
#ifdef DEBUG
    cout << "L is" << endl << L;
    cout << "invL.rows() =" << invL.rows() << " invL.cols() = " << invL.cols() << endl;
    cout << "H.rows() =" << H.rows() << " H.cols() = " << H.cols() << endl;
    cout << "invL_Ht.rows() =" << invL_Ht.rows() << " invL_Ht.cols() = " << invL_Ht.cols() << endl;
#endif        
        
        // invL_Ht = invL*H'
        Blas_Mat_Mat_Trans_Mult(invL,H,invL_Ht);
        
        // transformedObs = invL_Ht*observations
        Blas_Mat_Vec_Mult(invL_Ht,observations.col(iProcessedObservation),transformedObs);
        
        // we start by the root node
        iCurrentNode = 0;    
        
        while(nodes[iCurrentNode].height<_nInputs)
        {
            for(iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
            {
                // the parent node is replicated
                tTreeNode child = nodes[iCurrentNode];
                
                child.height++;
                child.symbolsVector(_nInputs-child.height) = _alphabet[iAlphabet];
                
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
            
#ifdef DEBUG
        cout << "symbols vector detected at " << iProcessedObservation << " = " << endl << nodes[iCurrentNode].symbolsVector;
#endif        
        
        // for the next iteration
        nodes.clear();
        
    } //for(int iProcessedObservation=_preambleLength;iProcessedObservation<_iLastSymbolVectorToBeDetected;iProcessedObservation++)
}

int KnownFlatChannelOptimalAlgorithm::iBestLeaf(const vector<tTreeNode> &nodes)
{
    int iBest = -1;
    double bestCost = 0.0;
    
#ifdef DEBUG3
    cout << "en iBestLeaf" << endl;
#endif
    
    for(uint i=0;i<nodes.size();i++)
    {
        // if it isn't a leaf node
        if(nodes[i].children.size()!=0)
            continue;
        
#ifdef DEBUG3
        cout << "iBest = " << iBest << " bestCost = " << bestCost << " nodes[i].cost = " << nodes[i].cost << endl;
#endif        
        
        if(iBest==-1 || nodes[i].cost < bestCost)
        {
            iBest = i;
            bestCost = nodes[i].cost;
        }
    }
    
    return iBest;
}


tMatrix KnownFlatChannelOptimalAlgorithm::getDetectedSymbolVectors()
{
    return _detectedSymbols;
}
