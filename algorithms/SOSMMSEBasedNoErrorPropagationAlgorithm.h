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
#ifndef SOSMMSEBASEDNOERRORPROPAGATIONALGORITHM_H
#define SOSMMSEBASEDNOERRORPROPAGATIONALGORITHM_H

#include <SOSMMSEBasedAlgorithm.h>

class SOSMMSEBasedNoErrorPropagationAlgorithm: public SOSMMSEBasedAlgorithm
{
protected:
	const MatrixXd _trueSymbols;
	const MIMOChannel * const _channel;
public:
    SOSMMSEBasedNoErrorPropagationAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, KalmanEstimator* kalmanEstimator, MatrixXd preamble, uint smoothingLag, SOSMMSEDetector *linearDetector, std::vector<double> ARcoefficients, const MatrixXd &symbols, const MIMOChannel * const channel):SOSMMSEBasedAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, kalmanEstimator, preamble, smoothingLag, linearDetector, ARcoefficients, true),_trueSymbols(symbols), _channel(channel)
	{
	}
	
	virtual MatrixXd getPreviousInterferingSymbols(uint iCurrentObservation) {return _trueSymbols.block(0,iCurrentObservation-_channelOrder+1,_nInputs,_channelOrder-1);}
	
	virtual std::vector<MatrixXd> getChannelMatricesForInterferenceCancellation(std::vector<MatrixXd> ARmatricesBuffer, uint iObservation) const
	{
		std::vector<MatrixXd> matrices(_d+1);
		
		for(uint i=0;i<=_d;i++)
			matrices[i] = _channel->at(iObservation+i);
		
		return matrices;
			
	}
};

#endif