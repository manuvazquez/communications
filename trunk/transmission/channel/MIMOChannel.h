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
    int nInputsnOutputsMemory(int n) const {return _nInputs*_nOutputs*memory(n);};
    int nInputsMemory(int n) const {return _nInputs*memory(n);}
	
	//! It returns the memory of the (possibly time-varying) channel
	/*!
		/param n time instant
		/return the memory of the channel at the given time instant
	*/
    virtual int memory(int n) const = 0;
	
	//! It returns the maximum memory of the channel
    virtual int effectiveMemory() const = 0;
    
    virtual MatrixXd at(int n) const = 0;
    
    virtual MatrixXd getTransmissionMatrix(const int n) const { return at(n);}
    
    MatrixXd transmit(const MatrixXd &symbols,const Noise &noise) const;
    
    vector<MatrixXd> range(int a,int b);
};

#endif
