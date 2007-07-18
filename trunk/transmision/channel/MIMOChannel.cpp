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

MIMOChannel::MIMOChannel(int nTx,int nRx,int length):_nTx(nTx),_nRx(nRx),_length(length),_nTxnRx(_nTx*_nRx)
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
	if(symbols.rows()!=_nTx)
		throw RuntimeException("MIMOChannel::Transmit: symbol vectors length is wrong.");
	else if(symbols.cols()>_length)
		throw RuntimeException("MIMOChannel::Transmit: channel length is shorter than then number of symbol vectors.");
	else if(noise.Nr()!=_nRx || symbols.cols()>noise.Length())
		throw RuntimeException("MIMOChannel::Transmit: missmatched noise dimensions.");

	// the number of resulting observations depends on the channel _memory
	int nObservations = symbols.cols() - (EffectiveMemory() - 1);

	if(nObservations<1)
		throw RuntimeException("MIMOChannel::Transmit: not enough symbol vectors for this channel _memory.");

	tMatrix observations = tMatrix(_nRx,symbols.cols());

	int j;
	tVector currentObservationVector(_nRx);

	tRange allChannelMatrixRows(0,_nRx-1);

	for(int iSymbolVector=EffectiveMemory()-1;iSymbolVector<symbols.cols();iSymbolVector++)
	{
		// just for the sake of clarity
		tMatrix &currentChannelMatrix = (*this)[iSymbolVector];

		//currentObservationVector will accumulate the contributions of the
		// different symbol vectors that participate in the current observation
		// (_memory >= 1)
		currentObservationVector = 0.0;

		for(j=0;j<Memory(iSymbolVector);j++)
		{
			// currentObservationVector = currentObservationVector + currentChannelMatrix(allChannelMatrixRows,*(new tRange(j*_nTx,(j+1)*_nTx-1)))*symbols.col(iSymbolVector-_memory+1+j)
			tRange rowsRange(j*_nTx,(j+1)*_nTx-1);
			Blas_Mat_Vec_Mult(currentChannelMatrix(allChannelMatrixRows,rowsRange),symbols.col(iSymbolVector-Memory(iSymbolVector)+1+j), currentObservationVector,1.0,1.0);
		}

		// the noise is added:
		//currentObservationVector = currentObservationVector + noise[iSymbolVector]
		Util::Add(currentObservationVector,noise[iSymbolVector],currentObservationVector);

		// the just computed observation is set in the observations matrix
		for(j=0;j<_nRx;j++)
			observations(j,iSymbolVector) = currentObservationVector(j);
	}

	return observations;
}

vector<tMatrix> MIMOChannel::Range(int a,int b)
{
    int nMatrices = b - a + 1;

    if(nMatrices<1)
        throw RuntimeException("MIMOChannel::Range: selected range of time is invalid.");
    vector<tMatrix> res(nMatrices);

    for(int i=0;i<nMatrices;i++)
        res[i] = operator[](a+i);

    return res;
}
