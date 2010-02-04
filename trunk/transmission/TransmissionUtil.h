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
#ifndef TRANSMISSIONUTIL_H
#define TRANSMISSIONUTIL_H

/**
    @author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <vector>
#include <Bits.h>
#include <Alphabet.h>
#include <Util.h>

class TransmissionUtil{
private:
    static void BERComputingChecks(const Bits &bits1,int from1,int to1,const Bits &bits2,int from2,int to2);
public:
    static double computeBER(const Bits &bits1,int from1,int to1,const Bits &bits2,int from2,int to2);
//     static double computeSER(const MatrixXd &sourceSymbols,const MatrixXd &detectedSymbols,const vector<vector<bool> > &mask,vector<vector<uint> > permutations,const Alphabet * const alphabet);
    static double computeBERsolvingAmbiguity(const Bits &sourceBits,int from1,int to1,const Bits &detectedBits,int from2,int to2,vector<vector<uint> > permutations);
    static VectorXd MSEalongTime(const std::vector<MatrixXd> &estimatedChannelMatrices,uint from1,uint to1,const std::vector<MatrixXd> &trueChannelMatrices,uint from2,uint to2);
    static MatrixXd generateTrainingSequence(const Alphabet &alphabet,uint nInputs,uint length);
};

#endif
