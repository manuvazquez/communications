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
