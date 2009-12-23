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
#ifndef BITS_H
#define BITS_H

/**
    @author Manu <manu@rustneversleeps>
*/

#include "types.h"
#include "exceptions.h"
#include <Random.h>
#include <vector>

class Bits{

private:
    uint _nStreams, _nBitsByStream,_nBits;
    tBit *_matrix;

public:
    Bits();
    Bits(int nStreams, int nBitsByStream,Random &randomGenerator);
    Bits(tBit *matrix,int nStreams,int nBitsByStream);
    Bits& operator=(const Bits& bits);
    Bits(const Bits& bits);
    ~Bits();

    void print() const;
    Bits differentialEncoding();
    Bits differentialDecoding();
    tBit operator()(int i,int j) const {return _matrix[i*_nBitsByStream+j];}
    int nStreams() const { return _nStreams;}
    int nBitsPerStream() const {return _nBitsByStream;}
    bool operator==(const Bits &bits) const;
    std::vector<tBit> GetStream(int index) const;
    void inject(int index,const std::vector<tBit> &stream);
    void invertStream(int index);

    // returns the number of non coincident bits
    int operator-(const Bits &bits) const;
};

#endif
