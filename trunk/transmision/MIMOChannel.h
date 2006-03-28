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
#include <excepcionesTransmision.h>
#include <Util.h>

using namespace la;

class MIMOChannel{
protected:
	int nTx, nRx, memory,length,nTx_nRx,nTx_nRx_memory,nTx_memory;

public:
    MIMOChannel();
	MIMOChannel(int nTx,int nRx, int memory, int length);
    ~MIMOChannel();

	int Nt() { return nTx;};
	int Nr() { return nRx;};
	int Memory() {return memory;};
	int Length() {return length;};
	int NtNr() {return nTx_nRx;};
	int NtNrMemory() {return nTx_nRx_memory;};
	int NtMemory() {return nTx_memory;};
	virtual tMatrix& operator[](int n) = 0;
	tMatrix Transmit(tMatrix &symbols,Noise &noise);
};

#endif
