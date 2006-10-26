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
#include <iostream>
#include <time.h>
#include "Bits.h"
// #include <Random.h>

Bits::Bits()
{
	nStreams = 0;
	nBitsByStream = 0;
	matrix = NULL;
}

Bits::Bits(int nStreams, int nBitsByStream,Random &randomGenerator)
{
	this->nStreams = nStreams;
	this->nBitsByStream = nBitsByStream;

// 	Random randomGenerator(1234);

	matrix = new tBit[nStreams*nBitsByStream];
	for(int i=nStreams*nBitsByStream;i--;)
	{
		matrix[i] = randomGenerator.randn() > 0 ? 1 : 0;
	}
}

Bits::Bits(tBit *matrix,int nStreams,int nBitsByStream): nStreams(nStreams),nBitsByStream(nBitsByStream),matrix(matrix)
{
}

Bits& Bits::operator=(const Bits& bits)
{
	cout << "Calling = operator..." << endl;
	nStreams = bits.nStreams;
	nBitsByStream = bits.nBitsByStream;
	if(matrix!=NULL)
	{
		delete[] matrix;
// 		cout << "Not null!!" << endl;
	}
	matrix = new tBit[bits.nStreams*bits.nBitsByStream];
	for(int i=nStreams*nBitsByStream;i--;)
		matrix[i] = bits.matrix[i];
	return *this;
}

Bits::Bits(const Bits& bits):
nStreams(bits.nStreams),nBitsByStream(bits.nBitsByStream)
{
	cout << "Calling copy constructor..." << endl;
	matrix = new tBit[bits.nStreams*bits.nBitsByStream];
	for(int i=nStreams*nBitsByStream;i--;)
		matrix[i] = bits.matrix[i];
}

Bits::~Bits()
{
	delete[] matrix;
}

void Bits::Print()
{
	int i,j;

	for(i=0;i<nStreams;i++)
	{
		for(j=0;j<nBitsByStream;j++)
			cout << matrix[i*nBitsByStream+j];
		cout << endl;
	}
}

Bits Bits::DifferentialEncoding()
{
	//Bits res(nStreams,nBitsByStream+1);
	Bits res;
	res.nStreams = nStreams;
	res.nBitsByStream = nBitsByStream+1;
	res.matrix = new tBit[res.nStreams*res.nBitsByStream];

	int i,j;
	for(i=0;i<nStreams;i++)
	{
		// differential encoding assumes the first bits equals 0
		res.matrix[i*res.nBitsByStream] = 0;
		for(j=0;j<nBitsByStream;j++)
			res.matrix[i*res.nBitsByStream+j+1] = (res.matrix[i*res.nBitsByStream+j] + matrix[i*nBitsByStream+j]) % 2;
	}
	return res;
}

Bits Bits::DifferentialDecoding()
{
	if(nBitsByStream<2)
		throw RuntimeException("2 bits by stream needed at least for differential decoding.");

	Bits res;
	res.nStreams = nStreams;
	res.nBitsByStream = nBitsByStream-1;
	res.matrix = new tBit[res.nStreams*res.nBitsByStream];

	int i,j;
	for(i=0;i<nStreams;i++)
	{
		for(j=1;j<nBitsByStream;j++)
		{
// 			cout << "i " << i << " j " << j << " res.nBitsByStream=" << res.nBitsByStream << endl;
			res.matrix[i*res.nBitsByStream+j-1] = (matrix[i*nBitsByStream+j-1] + matrix[i*nBitsByStream+j]) % 2;
		}
	}
	return res;
}

bool Bits::operator==(const Bits &bits) const
{
	for(int i=nStreams*nBitsByStream;i--;)
		if(matrix[i]!=bits.matrix[i])
			return false;
	return true;
}

int Bits::operator-(const Bits &bits) const
{
	int res = 0;
	for(int i=nStreams*nBitsByStream;i--;)
		if(matrix[i]!=bits.matrix[i])
			res++;
	return res;
}
