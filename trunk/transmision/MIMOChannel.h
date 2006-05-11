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
#ifndef MIMOCHANNEL_H
#define MIMOCHANNEL_H

/**
	@author Manu <manu@rustneversleeps>
*/

// #include <lapackpp/gmd.h>
#include <types.h>
#include <Noise.h>
#include <exceptions.h>
#include <Util.h>

using namespace la;

class MIMOChannel{
protected:
	int _nTx, _nRx, _memory,_length,_nTxnRx,_nTxnRxMemory,_nTxMemory;

public:
    MIMOChannel();
	MIMOChannel(int nTx,int nRx, int memory, int length);

	int Nt() { return _nTx;};
	int Nr() { return _nRx;};
	int Memory() {return _memory;};
	int Length() {return _length;};
	int NtNr() {return _nTxnRx;};
	int NtNrMemory() {return _nTxnRxMemory;};
	int NtMemory() {return _nTxMemory;};
	virtual tMatrix& operator[](int n) = 0;
	tMatrix Transmit(tMatrix &symbols,Noise &noise);
};

#endif
