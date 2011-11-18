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

// so that eigen doesn't complain about A!=B matrix operations (it's been told that it's a bug...)
#include<Eigen/Core>

#include<defines.h>

MIMOChannel::MIMOChannel(uint nInputs,uint nOutputs,uint length):_nInputs(nInputs),_nOutputs(nOutputs),_length(length),_nInputsnOutputs(_nInputs*_nOutputs)
{
}

/**
 * Transmits...
 * @param symbols
 * @param noise
 * @return
 */
 MatrixXd MIMOChannel::transmit(const MatrixXd &symbols,const Noise &noise) const
{
    if(symbols.rows()!=_nInputs)
        throw RuntimeException("MIMOChannel::transmit: symbol vectors length is wrong.");
    else if(symbols.cols()>_length)
        throw RuntimeException("MIMOChannel::transmit: channel length is shorter than then number of symbol vectors.");
    else if(noise.nOutputs()!=_nOutputs || symbols.cols()>noise.length())
        throw RuntimeException("MIMOChannel::transmit: missmatched noise dimensions.");

    // the number of resulting observations depends on the channel _memory
    int nObservations = symbols.cols() - (effectiveMemory() - 1);

    if(nObservations<1)
        throw RuntimeException("MIMOChannel::transmit: not enough symbol vectors for this channel _memory.");

	MatrixXd observations = FUNNY_VALUE * MatrixXd::Ones(_nOutputs,symbols.cols()).array();

    int j;
    
    for(int iSymbolVector=effectiveMemory()-1;iSymbolVector<symbols.cols();iSymbolVector++)
    {
        // just for the sake of clarity
        MatrixXd currentChannelMatrix = getTransmissionMatrix(iSymbolVector);

        //we will accumulate the contributions of the different symbol vectors that participate in the current observation (_memory >= 1).
        // Besides, we will always accumulate the noise.
        observations.col(iSymbolVector) = noise.at(iSymbolVector);

        for(j=0;j<memory(iSymbolVector);j++)
            observations.col(iSymbolVector) += currentChannelMatrix.block(0,j*_nInputs,_nOutputs,_nInputs)*symbols.col(iSymbolVector-memory(iSymbolVector)+1+j);
    }

    return observations;
}

vector<MatrixXd> MIMOChannel::range(int a,int b)
{
    int nMatrices = b - a + 1;

    if(nMatrices<1)
        throw RuntimeException("MIMOChannel::range: selected range of time is invalid.");

    vector<MatrixXd> res(nMatrices);

    for(int i=0;i<nMatrices;i++)
        res[i] = at(a+i);

    return res;
}

std::vector<uint> MIMOChannel::getInputsZeroCrossings(uint iFrom, uint length) const
{
  MatrixXd lastSignsMatrix = Util::sign(at(iFrom));
  
  std::vector<uint> res;
  res.push_back(iFrom);
  
  uint iChannelMatrix;
  
  for(iChannelMatrix=iFrom+1;iChannelMatrix<iFrom+length;iChannelMatrix++)
  {
	if(Util::sign(at(iChannelMatrix))!=lastSignsMatrix)
	{
	  lastSignsMatrix = Util::sign(at(iChannelMatrix));
	  res.push_back(iChannelMatrix);
	}
  }

  res.push_back(iChannelMatrix);

  return res;
}