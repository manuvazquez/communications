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

#include <string>
#include <types.h>
#include <Noise.h>
#include <exceptions.h>
#include <Util.h>

class MIMOChannel{
  
protected:
    uint _nInputs, _nOutputs,_length,_nInputsnOutputs;
public:
	MIMOChannel();
    MIMOChannel(uint nInputs,uint nOutputs,uint length);
    virtual ~MIMOChannel() {};
	
	virtual std::string name() const = 0;
    
    virtual uint nInputs() const { return _nInputs;}
    virtual uint nOutputs() const { return _nOutputs;}

	/*!
	  It returns the number of rows of the REAL channel matrix that represents the channel (it is usually the number of outputs)
	  \return number of rows of the internal channel coefficients matrix
	*/
	virtual uint channelCoefficientsMatrixRows() const { return _nOutputs;}

	/*!
	  It returns the number of columns of the REAL channel matrix that represents the channel (it is usually the number of inputs times the channel order)
	  \return number of columns of the internal channel coefficients matrix
	*/
	virtual uint channelCoefficientsMatrixCols() const { return _nInputs*effectiveMemory();}

	virtual uint length() const {return _length;};
    virtual uint nInputsnOutputs() const {return _nInputsnOutputs;};
    virtual uint nInputsnOutputsMemory(uint n) const {return _nInputs*_nOutputs*memory(n);};
    virtual uint nInputsMemory(uint n) const {return _nInputs*memory(n);}
	
	//! It returns the (possibly time-varying) memory of the channel
	/*!
		\param n time instant
		\return the memory of the channel at the given time instant
	*/
    virtual uint memory(uint n) const = 0;
	
	//! It returns the maximum memory of the channel
    virtual uint effectiveMemory() const = 0;

	//! It returns a matrix representing the channel at a given time instant
	/*!
		\param n time instant
		\return matrix that represents the channel at time \ref n
	*/
    virtual MatrixXd at(uint n) const = 0;
    
	//! It returns the channel matrix at a given time instant, as used in the signal model (it could be multiplied by a codes matrix in CDMA, e.g.)
	/*!
		\param n time instant
		\return channel matrix at time \ref n
	*/
    virtual MatrixXd getTransmissionMatrix(const uint n) const { return at(n);}
    
    virtual MatrixXd transmit(const MatrixXd &symbols,const Noise &noise) const;
    
    virtual std::vector<MatrixXd> range(uint a,uint b);
	
	std::vector<MatrixXd> getChannelMatrices();
	
	//! It sets the matrix representing the channel at the given time instant
	/*!
		\param n time instant
		\param mat matrix to be set
	*/
    virtual void set(uint n, MatrixXd mat)
    {
	  throw RuntimeException("MIMOChannel::set: the method is not implemented in the subclass.");
	}
};

#endif
