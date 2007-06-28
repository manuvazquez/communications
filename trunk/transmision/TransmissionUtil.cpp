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

// #define DEBUG

using namespace std;

void TransmissionUtil::BERComputingChecks(const Bits &bits1,int from1,int to1,const Bits &bits2,int from2,int to2)
{
    if((to1-from1)!=(to2-from2))
   	{
   		cout << "Range 1: " << (to1-from1) << " | " << "Range 2: " << to2-from2 << endl;
        throw RuntimeException("BERComputingChecks: comparisons range length are different.");
    }

    if(to1<from1)
        throw RuntimeException("BERComputingChecks: comparisons range are negatives.");

    if(to1>bits1.NbitsByStream() || to2>bits2.NbitsByStream() || from1<0 || from2<0)
    {
    	cout << "bits1.NbitsByStream(): " << bits1.NbitsByStream() << ",bits2.NbitsByStream(): " << bits2.NbitsByStream() << endl;
    	cout << "from1: " << from1 << ", to1: " << to1 << " | " << "from2: " << from2 << ", to2: " << to2 << endl;
        throw RuntimeException("BERComputingChecks: one or several comparison limits are wrong.");
    }

    if(bits1.Nstreams()!=bits2.Nstreams())
        throw RuntimeException("BERComputingChecks: bits objects have different number of streams.");
}

double TransmissionUtil::ComputeBER(const Bits &bits1,int from1,int to1,const Bits &bits2,int from2,int to2)
{
	if(bits2.NbitsByStream()==0 || bits2.Nstreams()==0)
		return 0.0;

	BERComputingChecks(bits1,from1,to1,bits2,from2,to2);

	int length = to1-from1;
    int errors = 0;

    for(int iBits1=from1,iBits2=from2;iBits1<to1;iBits1++,iBits2++)
        for(int iStream=0;iStream<bits1.Nstreams();iStream++)
            if(bits1(iStream,iBits1)!=bits2(iStream,iBits2))
                errors++;

    return (double)errors/(double)(length*bits1.Nstreams());
}

double TransmissionUtil::ComputeBERsolvingAmbiguity(const Bits &sourceBits,int from1,int to1,const Bits &detectedBits,int from2,int to2,vector<vector<uint> > permutations)
{
	if(detectedBits.NbitsByStream()==0 || detectedBits.Nstreams()==0)
		return 0.0;

	BERComputingChecks(sourceBits,from1,to1,detectedBits,from2,to2);

	int length = to1-from1;

	// max number of errors is length*sourceBits.Nstreams()
	int minErrors = length*sourceBits.Nstreams()+1;

    for(uint iPermut=0;iPermut<permutations.size();iPermut++)
    {
    	int errorsPermutation = 0;

        for(uint iStream=0;iStream<permutations[iPermut].size();iStream++)
        {
        	int errorsInverting=0,errorsWithoutInverting=0;

        	// without inverting bits
        	for(int iSourceStream=from1,iDetectedStream=from2;iSourceStream<to1;iSourceStream++,iDetectedStream++)
        		errorsWithoutInverting += (sourceBits(iStream,iSourceStream) != detectedBits(permutations[iPermut][iStream],iDetectedStream));

        	// inverting bits
        	for(int iSourceStream=from1,iDetectedStream=from2;iSourceStream<to1;iSourceStream++,iDetectedStream++)
        		errorsInverting += (sourceBits(iStream,iSourceStream) == detectedBits(permutations[iPermut][iStream],iDetectedStream));

        	if(errorsWithoutInverting<errorsInverting)
        		errorsPermutation += errorsWithoutInverting;
        	else
        		errorsPermutation += errorsInverting;
        }

        if(errorsPermutation<minErrors)
        	minErrors = errorsPermutation;
    }

    return (double)minErrors/(double)(length*sourceBits.Nstreams());
}

tVector TransmissionUtil::MSEalongTime(const std::vector<tMatrix> &estimatedChannelMatrices,int from1,int to1,const std::vector<tMatrix> &trueChannelMatrices,int from2,int to2)
{

	tVector res(to1-from1+1);

    // if the algorithm didn't make channel estimation
	if(estimatedChannelMatrices.size()==0)
	{
		res = 0.0;
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

#ifdef DEBUG
	cout << "antes de inicializar res" << endl;
#endif

#ifdef DEBUG
	cout << "hola" << endl;
#endif

	// if the channel is Sparkling memory, the channel matrices of the real channel may have different sizes
	try {
		for(int iSource1=from1,iSource2=from2,iRes=0;iSource1<=to1;iSource1++,iSource2++,iRes++)
		{
			// the square error committed by the estimated matrix is normalized by the squared Frobenius norm (i.e. the sum of all the elements squared) of the real channel matrix
			res(iRes) = Util::SquareErrorPaddingWithZeros(trueChannelMatrices.at(iSource2),estimatedChannelMatrices.at(iSource1))/pow(Blas_NormF(trueChannelMatrices.at(iSource2)),2.0);
#ifdef DEBUG
// 			cout << "comparando" << endl << trueChannelMatrices.at(iSource2) << "y" << endl << estimatedChannelMatrices.at(iSource1) << endl;
			cout << "res(" << iRes << ") = " << res(iRes) << " res.size() = " << res.size() << endl;
#endif
		}
	} catch (IncompatibleOperandsException) {
		return res;
	}
#ifdef DEBUG
	cout << "res.zie " << res.size() << endl;
#endif

	return res;
}
