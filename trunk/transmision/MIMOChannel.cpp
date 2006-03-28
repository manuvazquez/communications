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
#include "MIMOChannel.h"
#include <lapackpp/blas2pp.h>

using namespace la;

MIMOChannel::MIMOChannel()
{
	nTx = 0;
	nRx = 0;
	memory = 0;
	length = 0;
	nTx_nRx = 0;
	nTx_nRx_memory = 0;
	nTx_memory = 0;
}

MIMOChannel::MIMOChannel(int nTx,int nRx, int memory, int length)
{
	this->nTx = nTx;
	this->nRx = nRx;
	this->memory = memory;
	this->length = length;
	this->nTx_nRx = nTx*nRx;
	this->nTx_nRx_memory = nTx*nRx*memory;
	this->nTx_memory = nTx*memory;
}

MIMOChannel::~MIMOChannel()
{
}

/**
 * Transmits...
 * @param symbols 
 * @param noise 
 * @return 
 */
tMatrix MIMOChannel::Transmit(tMatrix &symbols,Noise &noise)
{
	if(symbols.rows()!=nTx)
		throw RuntimeException("Symbol vectors length is wrong.");
	else if(symbols.cols()>length)
		throw RuntimeException("Channel length is shorter than then number of symbol vectors.");
	else if(noise.Nr()!=nRx || symbols.cols()>noise.Length())
		throw RuntimeException("Missmatched noise dimensions.");
	
	// the number of resulting observations depends on the channel memory
	int nObservations = symbols.cols() - (memory - 1);

	if(nObservations<1)
		throw RuntimeException("Not enough symbol vectors for this channel memory.");

	tMatrix observations = tMatrix(nRx,symbols.cols());

	int j;
	tVector currentObservationVector(nRx);

	tRange allChannelMatrixRows(0,nRx-1);

	for(int iSymbolVector=memory-1;iSymbolVector<symbols.cols();iSymbolVector++)
	{
		// just for the sake of clarity
		tMatrix &currentChannelMatrix = (*this)[iSymbolVector];

		//currentObservationVector will accumulate the contributions of the
		// different symbol vectors that participate in the current observation
		// (memory >= 1)
		currentObservationVector = 0.0;
		
		for(j=0;j<memory;j++)
		{
			// currentObservationVector = currentObservationVector + currentChannelMatrix(allChannelMatrixRows,*(new tRange(j*nTx,(j+1)*nTx-1)))*symbols.col(iSymbolVector-memory+1+j)
			Blas_Mat_Vec_Mult(currentChannelMatrix(allChannelMatrixRows,*(new tRange(j*nTx,(j+1)*nTx-1))),symbols.col(iSymbolVector-memory+1+j), currentObservationVector,1.0,1.0);
		}

		// the noise is added:
		//currentObservationVector = currentObservationVector + noise[iSymbolVector]
		Util::Add(currentObservationVector,noise[iSymbolVector],currentObservationVector);

		// the just computed observation is set in the observations matrix	
		for(j=0;j<nRx;j++)
			observations(j,iSymbolVector) = currentObservationVector(j);
	}
	
	return observations;
}
