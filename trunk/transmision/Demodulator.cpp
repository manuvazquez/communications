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
#include "Demodulator.h"

Demodulator::Demodulator()
{
}


Demodulator::~Demodulator()
{
}

Bits Demodulator::Demodulate(const tMatrix &symbols,Alphabet alphabet)
{
	int nBitsByStream = symbols.cols()*alphabet.NbitsBySymbol();
	int nStreams = symbols.rows();
	tBit *matrix = new tBit[nStreams*nBitsByStream];

	int iBit,j,k;

	for(int i=0;i<nStreams;i++)
	{
		iBit = 0;
		for(j=0;j<symbols.cols();j++)
		{
			vector<tBit> bitsSequence = alphabet[(tSymbol)symbols(i,j)];
			for(k=0;k<alphabet.NbitsBySymbol();k++,iBit++)
				matrix[i*nBitsByStream+iBit] = bitsSequence[k];
		}
	}
	return Bits(matrix,nStreams,nBitsByStream);
	
}
