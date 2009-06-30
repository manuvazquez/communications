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

Bits::Bits()
{
	_nStreams = 0;
	_nBitsByStream = 0;
	_nBits = 0;
	_matrix = NULL;
}

Bits::Bits(int nStreams, int nBitsByStream,Random &randomGenerator):_nStreams(nStreams),_nBitsByStream(nBitsByStream),_nBits(nStreams*nBitsByStream)
{
	_matrix = new tBit[_nBits];
	for(int i=_nBits;i--;)
	{
		_matrix[i] = randomGenerator.randn() > 0 ? 1 : 0;
	}
}

Bits::Bits(tBit *matrix,int nStreams,int nBitsByStream): _nStreams(nStreams),_nBitsByStream(nBitsByStream),_nBits(nStreams*nBitsByStream),_matrix(matrix)
{
}

Bits& Bits::operator=(const Bits& bits)
{
	if(_nStreams!=bits._nStreams || _nBitsByStream!=bits._nBitsByStream)
	{
		_nStreams = bits._nStreams;
		_nBitsByStream = bits._nBitsByStream;
		_nBits = bits._nBits;
		delete[] _matrix;
		_matrix = new tBit[_nBits];
	}

	for(int i=_nBits;i--;)
		_matrix[i] = bits._matrix[i];

	return *this;
}

Bits::Bits(const Bits& bits):
_nStreams(bits._nStreams),_nBitsByStream(bits._nBitsByStream),_nBits(bits._nBits)
{
	_matrix = new tBit[_nBits];
	for(int i=_nBits;i--;)
		_matrix[i] = bits._matrix[i];
}

Bits::~Bits()
{
	delete[] _matrix;
}

void Bits::Print() const
{
	uint i,j;

	for(i=0;i<_nStreams;i++)
	{
		for(j=0;j<_nBitsByStream;j++)
			cout << _matrix[i*_nBitsByStream+j];
		cout << endl;
	}
}

Bits Bits::DifferentialEncoding()
{
	Bits res;
	res._nStreams = _nStreams;
	res._nBitsByStream = _nBitsByStream+1;
	res._nBits = res._nStreams*res._nBitsByStream;
	res._matrix = new tBit[res._nBits];

	uint i,j;
	for(i=0;i<_nStreams;i++)
	{
		// differential encoding assumes the first bits equals 0
		res._matrix[i*res._nBitsByStream] = 0;
		for(j=0;j<_nBitsByStream;j++)
			res._matrix[i*res._nBitsByStream+j+1] = (res._matrix[i*res._nBitsByStream+j] + _matrix[i*_nBitsByStream+j]) % 2;
	}
	return res;
}

Bits Bits::DifferentialDecoding()
{
	if(_nBitsByStream<2)
		throw RuntimeException("2 bits by stream needed at least for differential decoding.");

	Bits res;
	res._nStreams = _nStreams;
	res._nBitsByStream = _nBitsByStream-1;
	res._nBits = res._nStreams*res._nBitsByStream;
	res._matrix = new tBit[res._nBits];

	uint i,j;
	for(i=0;i<_nStreams;i++)
	{
		for(j=1;j<_nBitsByStream;j++)
			res._matrix[i*res._nBitsByStream+j-1] = (_matrix[i*_nBitsByStream+j-1] + _matrix[i*_nBitsByStream+j]) % 2;
	}
	return res;
}

bool Bits::operator==(const Bits &bits) const
{
	for(int i=_nBits;i--;)
		if(_matrix[i]!=bits._matrix[i])
			return false;
	return true;
}

int Bits::operator-(const Bits &bits) const
{
	int res = 0;
	for(int i=_nStreams*_nBitsByStream;i--;)
		if(_matrix[i]!=bits._matrix[i])
			res++;
	return res;
}

vector<tBit> Bits::GetStream(int index) const
{
	vector<tBit> res(_nBitsByStream);

	for(uint i=index*_nBitsByStream,j=0;j<_nBitsByStream;i++,j++)
		res[j] = _matrix[i];

	return res;
}

void Bits::Inject(int index,const std::vector<tBit> &stream)
{
	if(stream.size()!=_nBitsByStream)
		throw RuntimeException("Bits::Inject: the stream has not the correct number of bits.");

	for(uint i=index*_nBitsByStream,j=0;j<_nBitsByStream;i++,j++)
		_matrix[i] = stream[j];
}

void Bits::InvertStream(int index)
{
	for(uint i=index*_nBitsByStream,j=0;j<_nBitsByStream;i++,j++)
		_matrix[i] = _matrix[i]^1;
}
