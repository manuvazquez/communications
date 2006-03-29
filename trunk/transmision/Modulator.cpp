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
#include "Modulator.h"

Modulator::Modulator()
{
}


Modulator::~Modulator()
{
}

tMatrix Modulator::Modulate(const Bits &bits, Alfabeto alphabet)
{
	if((bits.NbitsByStream()% alphabet.NbitsPorSimbolo())!=0)
		cout << "Too many bits." << endl;
	int nSymbolsByStream = bits.NbitsByStream()/ alphabet.NbitsPorSimbolo();

	tMatrix res(bits.Nstreams(),nSymbolsByStream);

	// once filled, it will converted to a symbol by alphabet
	vector<tBit> currentBitSequence(alphabet.NbitsPorSimbolo());

	int processedBits,j,iSymbol;
	for(int i=0;i<bits.Nstreams();i++)
	{
		processedBits = 0;
		iSymbol = 0;
		while((processedBits+alphabet.NbitsPorSimbolo()-1)<bits.NbitsByStream())
		{
			for(j=0;j<alphabet.NbitsPorSimbolo();j++)
				currentBitSequence[j] = bits(i,processedBits+j);

			processedBits += alphabet.NbitsPorSimbolo();

			res(i,iSymbol++) = (double) alphabet[currentBitSequence];
		}
	}
	return res;
}
