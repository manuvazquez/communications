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
#include "ChannelMatrixEstimator.h"

ChannelMatrixEstimator::ChannelMatrixEstimator(MatrixXd initialEstimation,uint N):_nOutputs(initialEstimation.rows()),_nChannelMatrixRows(initialEstimation.rows()),_nInputsXchannelOrder(initialEstimation.cols()),_nInputs(N),_nChannelCoeffs(initialEstimation.rows()*initialEstimation.cols()),_lastEstimatedChannelCoefficientsMatrix(initialEstimation)
{
    if(_nInputsXchannelOrder < _nInputs)
        throw RuntimeException("ChannelMatrixEstimator::ChannelMatrixEstimator: number of columns of \"initialEstimation\"  is less than N");

    // check erased because of "OneChannelOrderPerTransmitAtennaWrapperEstimator"
    /*
    if((_nInputsXchannelOrder % _nInputs) != 0)
        throw RuntimeException("ChannelMatrixEstimator::ChannelMatrixEstimator: number of columns of \"initialEstimation\"  is not a multiple of N");
    */

    if((_nInputsXchannelOrder % _nInputs) == 0)
        _channelOrder = _nInputsXchannelOrder/_nInputs;
    // _channelOrder=-1 accounts for the case of a "OneChannelOrderPerTransmitAtennaWrapperEstimator" being used, whose internal ChannelMatrixEstimator does not need to have a number of columns multiple of _nInputs
    else
        _channelOrder = 0;
}

uint ChannelMatrixEstimator::memory() const
{
    if(_channelOrder!=0)
        return _channelOrder;
    else
        throw RuntimeException("ChannelMatrixEstimator::Memory: this may not be a real channel matrix estimator: its number of columns is not a multiple of the number of transmitting antennas.");
}

vector<MatrixXd> ChannelMatrixEstimator::nextMatricesFromObservationsSequence(const MatrixXd &observations,vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo)
{
    if(observations.cols()<iTo)
        throw RuntimeException("ChannelMatrixEstimator::nextMatricesFromObservationsSequence: insufficient number of observations.");

    vector<MatrixXd> estimatedMatrices(iTo-iFrom);

    for(uint i=iFrom;i<iTo;i++)
        estimatedMatrices[i-iFrom] = nextMatrix(observations.col(i),symbolVectors.block(0,i-_channelOrder+1,_nInputs,_channelOrder),noiseVariances[i]);
    
    return estimatedMatrices;
}

std::vector<MatrixXd> ChannelMatrixEstimator::nextMatricesFromObservationsSequence(const MatrixXd &observations,std::vector<double> &noiseVariances,const MatrixXd &symbolVectors,uint iFrom,uint iTo,std::vector<MatrixXd> &channelEstimatesVariances)
{
	assert(observations.cols()>=iTo);
	assert(computesVariances());
	
    std::vector<MatrixXd> estimatedMatrices(iTo-iFrom);
	channelEstimatesVariances = std::vector<MatrixXd>(iTo-iFrom);

    for(uint i=iFrom;i<iTo;i++)
	{
        estimatedMatrices[i-iFrom] = nextMatrix(observations.col(i),symbolVectors.block(0,i-_channelOrder+1,_nInputs,_channelOrder),noiseVariances[i]);
		channelEstimatesVariances[i-iFrom] = getVariances();
	}
    
    return estimatedMatrices;
}

