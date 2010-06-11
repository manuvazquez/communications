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
#include "TransmissionUtil.h"

// #define PRINT_INFO

void TransmissionUtil::BERComputingChecks(const Bits &bits1,int from1,int to1,const Bits &bits2,int from2,int to2)
{
    if((to1-from1)!=(to2-from2))
    {
        cout << "Range 1: " << (to1-from1) << " | " << "Range 2: " << to2-from2 << endl;
        throw RuntimeException("BERComputingChecks: comparisons range length are different.");
    }

    if(to1<from1)
        throw RuntimeException("BERComputingChecks: comparisons ranges are negatives.");

    if(to1>bits1.nBitsPerStream() || to2>bits2.nBitsPerStream() || from1<0 || from2<0)
    {
        cout << "bits1.nBitsPerStream(): " << bits1.nBitsPerStream() << ",bits2.nBitsPerStream(): " << bits2.nBitsPerStream() << endl;
        cout << "from1: " << from1 << ", to1: " << to1 << " | " << "from2: " << from2 << ", to2: " << to2 << endl;
        throw RuntimeException("BERComputingChecks: one or several comparison limits are wrong.");
    }

    if(bits1.nStreams()!=bits2.nStreams())
        throw RuntimeException("BERComputingChecks: bits objects have different number of streams.");
}

double TransmissionUtil::computeBER(const Bits &bits1,int from1,int to1,const Bits &bits2,int from2,int to2)
{
    if(bits2.nBitsPerStream()==0 || bits2.nStreams()==0)
        return 0.0;

    BERComputingChecks(bits1,from1,to1,bits2,from2,to2);

    int length = to1-from1;
    int errors = 0;

    for(int iBits1=from1,iBits2=from2;iBits1<to1;iBits1++,iBits2++)
        for(int iStream=0;iStream<bits1.nStreams();iStream++)
            if(bits1(iStream,iBits1)!=bits2(iStream,iBits2))
                errors++;

    return (double)errors/(double)(length*bits1.nStreams());
}

double TransmissionUtil::computeBERsolvingAmbiguity(const Bits &sourceBits,int from1,int to1,const Bits &detectedBits,int from2,int to2,vector<vector<uint> > permutations)
{
    if(detectedBits.nBitsPerStream()==0 || detectedBits.nStreams()==0)
        return 0.0;

    BERComputingChecks(sourceBits,from1,to1,detectedBits,from2,to2);

    int length = to1-from1;

    uint iBestPermutation = 0;
    vector<int> bestPermutationSigns(sourceBits.nStreams(),1);

    // max number of errors is length*sourceBits.nStreams()
    int minErrors = length*sourceBits.nStreams()+1;

    bool bitsDiffer;
    
    for(uint iPermut=0;iPermut<permutations.size();iPermut++)
    {
        int errorsPermutation = 0;

        for(uint iStream=0;iStream<permutations[iPermut].size();iStream++)
        {
            int errorsInverting=0,errorsWithoutInverting=0;

            for(int iSourceStream=from1,iDetectedStream=from2;iSourceStream<to1;iSourceStream++,iDetectedStream++)
            {
                // do bits differ?
                bitsDiffer = (sourceBits(iStream,iSourceStream) != detectedBits(permutations[iPermut][iStream],iDetectedStream));
                
                // if they do, this entails an error
                errorsWithoutInverting += bitsDiffer;
                
                // or no error if the bit needs to be inverted due to the ambiguity
                errorsInverting += !bitsDiffer;
            }

            if(errorsWithoutInverting<errorsInverting)
            {
                errorsPermutation += errorsWithoutInverting;
                bestPermutationSigns[iStream] = 1;
            }
            else
            {
                errorsPermutation += errorsInverting;
                bestPermutationSigns[iStream] = -1;
            }
        }

        if(errorsPermutation<minErrors)
        {
            minErrors = errorsPermutation;
            iBestPermutation = iPermut;
        }
    }

#ifdef DEBUG
    cout << "dividing " << minErrors << " by " << length*sourceBits.nStreams() << endl;
    cout << "iBestPermutation = " << iBestPermutation << " and the signs" << endl;
    Util::print(bestPermutationSigns);
    Util::print(permutations);    
#endif
    return (double)minErrors/(double)(length*sourceBits.nStreams());
}

VectorXd TransmissionUtil::MSEalongTime(const std::vector<MatrixXd> &estimatedChannelMatrices,uint from1,uint to1,const std::vector<MatrixXd> &trueChannelMatrices,uint from2,uint to2)
{

    VectorXd res(to1-from1+1);

    // if the algorithm didn't make channel estimation
    if(estimatedChannelMatrices.size()==0)
    {
        res.setConstant(0.0);
        return res;
    }

    if((to1-from1)!=(to2-from2))
    {
        cout << "Range 1: " << (to1-from1) << " | " << "Range 2: " << to2-from2 << endl;
        throw RuntimeException("TransmissionUtil::MSEalongTime: comparisons range length are different.");
    }

    if(to1<from1)
        throw RuntimeException("TransmissionUtil::MSEalongTime: comparisons range are negatives.");

    if(to1>estimatedChannelMatrices.size() || to2>trueChannelMatrices.size() || from1<0 || from2<0)
    {
        cout << from1 << endl << to1 << endl << from2 << endl << to2 << endl << estimatedChannelMatrices.size() << endl << trueChannelMatrices.size() << endl;
        throw RuntimeException("TransmissionUtil::MSEalongTime: one or several comparison limits are wrong.");
    }

    // if the channel is Sparkling memory, the channel matrices of the real channel may have different sizes
    try {
        for(uint iSource1=from1,iSource2=from2,iRes=0;iSource1<=to1;iSource1++,iSource2++,iRes++)
        {
            // the square error committed by the estimated matrix is normalized by the squared Frobenius norm (i.e. the sum of all the elements squared) of the real channel matrix
            res(iRes) = Util::squareErrorPaddingWithZeros(trueChannelMatrices.at(iSource2),estimatedChannelMatrices.at(iSource1))/pow(trueChannelMatrices.at(iSource2).norm(),2.0);

        }
    } catch (IncompatibleOperandsException) {
        return res;
    }

    return res;
}

MatrixXd TransmissionUtil::generateTrainingSequence(const Alphabet &alphabet,uint nInputs,uint length)
{
    MatrixXd res = MatrixXd::Zero(nInputs,length);

    if((length % nInputs) != 0)
        throw RuntimeException("TransmissionUtil::GenerateTrainingSequence: length is not a multiple of the number of transmitting antennas.");

    for(uint i=0;i<length/nInputs;i++)
    {
        res(0,i*nInputs) = 1;
        res(1,i*nInputs) = 1;
        res(0,i*nInputs+1) = -1;
        res(1,i*nInputs+1) = 1;
    }
    return res;
}