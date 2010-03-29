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

Bits Demodulator::demodulate(const MatrixXd &symbols,Alphabet alphabet)
{
    int nBitsByStream = symbols.cols()*alphabet.nBitsPerSymbol();
    int nStreams = symbols.rows();
    tBit *matrix = new tBit[nStreams*nBitsByStream];

    int iBit,j,k;

    for(int i=0;i<nStreams;i++)
    {
        iBit = 0;
        for(j=0;j<symbols.cols();j++)
        {
            vector<tBit> bitsSequence = alphabet[(tSymbol)symbols(i,j)];
            for(k=0;k<alphabet.nBitsPerSymbol();k++,iBit++)
                matrix[i*nBitsByStream+iBit] = bitsSequence[k];
        }
    }
    return Bits(matrix,nStreams,nBitsByStream);
}

std::vector<std::vector<bool> > Demodulator::demodulate(const std::vector<std::vector<bool> > &mask,const Alphabet &alphabet)
{
  if(mask.size()==0)
	throw RuntimeException("Demodulator::demodulate: the mask is empty.");

  int nBitsPerStream = mask[0].size()*alphabet.nBitsPerSymbol();
  int nStreams = mask.size();
  std::vector<std::vector<bool> > res(nStreams,std::vector<bool>(nBitsPerStream));

  
  int iBit,k;
  uint j;

  for(int i=0;i<nStreams;i++)
  {
	iBit = 0;
	for(j=0;j<mask[i].size();j++)
	{
	  for(k=0;k<alphabet.nBitsPerSymbol();k++,iBit++)
		  res[i][iBit] = mask[i][j];
	}
  }
  return res;
}

Bits Demodulator::demodulate(const MatrixXd &symbols,const Alphabet &alphabet,const std::vector<std::vector<bool> > &mask)
{
//   cout << "la mascara es" << endl;
//   Util::print(mask);
//   cout << endl;
// 
//   cout << " y los simbolos" << endl << symbols << endl;
  
  if(mask.size()==0)
	throw RuntimeException("Demodulator::demodulate: the mask is empty.");
  
  if(symbols.rows()!=mask.size() || symbols.cols()!=mask[0].size())
	throw RuntimeException("Demodulator::demodulate: symbols and mask size doesn't match.");
  
  int nBitsPerStream = symbols.cols()*alphabet.nBitsPerSymbol();
  int nStreams = symbols.rows();
  tBit *matrix = new tBit[nStreams*nBitsPerStream];

  vector<tBit> noSymbolBitsSequence(alphabet.nBitsPerSymbol(),Bits::noDataBitValue());

  int iBit,j,k;

  for(int i=0;i<nStreams;i++)
  {
	iBit = 0;
	for(j=0;j<symbols.cols();j++)
	{
// 	  if(mask[i][j] && alphabet.doesItBelong(symbols(i,j)))
	  if(alphabet.doesItBelong(symbols(i,j)))
	  {
		vector<tBit> bitsSequence = alphabet[(tSymbol)symbols(i,j)];
		for(k=0;k<alphabet.nBitsPerSymbol();k++,iBit++)
			matrix[i*nBitsPerStream+iBit] = bitsSequence[k];
	  }else
		for(k=0;k<alphabet.nBitsPerSymbol();k++,iBit++)
			matrix[i*nBitsPerStream+iBit] = noSymbolBitsSequence[k];
	}
  }
  return Bits(matrix,nStreams,nBitsPerStream);
}