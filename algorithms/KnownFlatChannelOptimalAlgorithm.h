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
#ifndef KNOWNFLATCHANNELOPTIMALALGORITHM_H
#define KNOWNFLATCHANNELOPTIMALALGORITHM_H

#include <KnownChannelAlgorithm.h>

/**
It implements the optimal algorithm to detect in a flat channel by means of a tree search.

    @author Manu <manu@rustneversleeps>
*/

#include <Eigen/LU> 
#include <Eigen/Cholesky>

class KnownFlatChannelOptimalAlgorithm : public KnownChannelAlgorithm
{
private:
    typedef struct{
        double cost;
        std::vector<uint> children;
        uint height,id;
        VectorXd symbolsVector;
    } tTreeNode;
    
    uint iBestLeaf(const std::vector< KnownFlatChannelOptimalAlgorithm::tTreeNode >& nodes);
protected:
    const uint _preambleLength;
    MatrixXd _detectedSymbols;
    Alphabet *_extendedAlphabet;
    
    virtual const Alphabet *getAlphabetAt(uint time, int leafHeight) const { return _extendedAlphabet;}
public:
    KnownFlatChannelOptimalAlgorithm(string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, const MIMOChannel& channel, uint preambleLength);

    ~KnownFlatChannelOptimalAlgorithm();

    void run(MatrixXd observations, vector< double > noiseVariances);
    MatrixXd getDetectedSymbolVectors();

};

#endif
