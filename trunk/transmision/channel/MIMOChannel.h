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
    int _nInputs, _nOutputs,_length,_nInputsnOutputs;
public:
    MIMOChannel(int nInputs,int nOutputs,int length);
    virtual ~MIMOChannel() {};
    
    int nInputs() const { return _nInputs;};
    int nOutputs() const { return _nOutputs;};
    int length() const {return _length;};
    int nInputsnOutputs() const {return _nInputsnOutputs;};
    int nInputsnOutputsMemory(int n) const {return _nInputs*_nOutputs*Memory(n);};
    int nInputsMemory(int n) const {return _nInputs*Memory(n);};
    virtual int Memory(int n) const = 0;
    virtual int Effectivememory() const = 0;
    virtual tMatrix operator[](int n) const = 0;
    tMatrix transmit(tMatrix &symbols,Noise &noise);
    vector<tMatrix> range(int a,int b);
};

#endif
