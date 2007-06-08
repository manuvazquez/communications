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

#include <types.h>
#include <Noise.h>
#include <exceptions.h>
#include <Util.h>

class MIMOChannel{
protected:
	int _nTx, _nRx,_length,_nTxnRx;
public:
	MIMOChannel(int nTx,int nRx,int length);
	virtual ~MIMOChannel() {};

	int Nt() const { return _nTx;};
	int Nr() const { return _nRx;};
	int Length() const {return _length;};
	int NtNr() const {return _nTxnRx;};
	int NtNrMemory(int n) const {return _nTx*_nRx*Memory(n);};
	int NtMemory(int n) const {return _nTx*Memory(n);};
	virtual int Memory(int n) const = 0;
	virtual int EffectiveMemory() const = 0;
	virtual tMatrix& operator[](int n) const = 0;
	tMatrix Transmit(tMatrix &symbols,Noise &noise);
    vector<tMatrix> Range(int a,int b);
};

#endif
