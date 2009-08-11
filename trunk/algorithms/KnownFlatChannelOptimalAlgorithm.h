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
class KnownFlatChannelOptimalAlgorithm : public KnownChannelAlgorithm
{
private:
    typedef struct{
        double cost;
        std::vector<int> children;
        uint height,id;
        tVector symbolsVector;
    } tTreeNode;
    
    int iBestLeaf(const std::vector<tTreeNode> &nodes);
protected:
    const int _preambleLength;
    tMatrix _detectedSymbols;
    Alphabet *_extendedAlphabet;
    
    virtual const Alphabet *getAlphabetAt(int time, int leafHeight) const { return _extendedAlphabet;}
public:
    KnownFlatChannelOptimalAlgorithm(string name, Alphabet alphabet, int L, int Nr, int N, int iLastSymbolVectorToBeDetected, const MIMOChannel& channel, int preambleLength);

    ~KnownFlatChannelOptimalAlgorithm();

    void run(tMatrix observations, vector< double > noiseVariances);
    tMatrix getDetectedSymbolVectors();

};

#endif