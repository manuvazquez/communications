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

Bits::Bits(int nStreams, int nBitsByStream)
{
	this->nStreams = nStreams;
	this->nBitsByStream = nBitsByStream;
	
	Random randomGenerator(1234);

	matrix = new tBit[nStreams*nBitsByStream];
	int aux;
	float numGenerado;
	for(int i=nStreams*nBitsByStream;i--;)
		matrix[i] = randomGenerator.randn() > 0 ? 1 : 0;
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

Bits::~Bits()
{
}


