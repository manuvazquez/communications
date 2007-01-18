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
#include "OneChannelOrderPerTransmitAtennaWrapperEstimator.h"

// #define DEBUG

using namespace std;

OneChannelOrderPerTransmitAtennaWrapperEstimator::OneChannelOrderPerTransmitAtennaWrapperEstimator(tMatrix initialEstimation, int N, const vector<int> &antennasChannelOrders,ChannelMatrixEstimator *realEstimator): ChannelMatrixEstimator(initialEstimation, N),_antennasChannelOrders(antennasChannelOrders),_realEstimator(realEstimator),_involvedSymbolsVector(realEstimator->Cols())
{
	#ifdef DEBUG
		cout << "_antennasChannelOrders.size(): " << _antennasChannelOrders.size() << " _N: " << _N << endl;
		cout << "El tamaño de involvedSymb... es " << _involvedSymbolsVector.size() << endl;
	#endif
	if(_antennasChannelOrders.size()!=_N)
		throw RuntimeException("OneChannelOrderPerTransmitAtennaWrapperEstimator::OneChannelOrderPerTransmitAtennaWrapperEstimator: antennasChannelOrders size is wrong.");

	if(_lastEstimatedChannelMatrix.cols() != (_antennasChannelOrders[Util::Max(_antennasChannelOrders)]*_N))
		throw RuntimeException("OneChannelOrderPerTransmitAtennaWrapperEstimator::OneChannelOrderPerTransmitAtennaWrapperEstimator: the length of the antennas channel orders vector is not coherent with the initial estimation channel matrix dimensions");

	if(_realEstimator->Cols()!=Util::Sum(_antennasChannelOrders))
		throw RuntimeException("OneChannelOrderPerTransmitAtennaWrapperEstimator::OneChannelOrderPerTransmitAtennaWrapperEstimator: the underlying estimator is not coherent with \"antennasChannelOrders\".");
}

OneChannelOrderPerTransmitAtennaWrapperEstimator::OneChannelOrderPerTransmitAtennaWrapperEstimator(const OneChannelOrderPerTransmitAtennaWrapperEstimator &estimator):ChannelMatrixEstimator(estimator),_antennasChannelOrders(estimator._antennasChannelOrders),_involvedSymbolsVector(estimator._involvedSymbolsVector),_realEstimator(estimator._realEstimator->Clone())
{
}

OneChannelOrderPerTransmitAtennaWrapperEstimator::~OneChannelOrderPerTransmitAtennaWrapperEstimator()
{
	delete _realEstimator;
}


ChannelMatrixEstimator* OneChannelOrderPerTransmitAtennaWrapperEstimator::Clone()
{
	return new OneChannelOrderPerTransmitAtennaWrapperEstimator(*this);
}

double OneChannelOrderPerTransmitAtennaWrapperEstimator::Likelihood(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance)
{
	tVector symbolsVector = Util::ToVector(symbolsMatrix,columnwise);

	// symbols not involved in the observation are supressed in "_involvedSymbolsVector"
	OneChannelOrderPerTransmitAtennaMIMOChannel::CompleteSymbolsVectorToOnlyInvolvedSymbolsVector(symbolsVector,_N,_antennasChannelOrders,_involvedSymbolsVector);

	return _realEstimator->Likelihood(observations,_involvedSymbolsVector,noiseVariance);
}

tMatrix OneChannelOrderPerTransmitAtennaWrapperEstimator::NextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
	tVector symbolsVector = Util::ToVector(symbolsMatrix,columnwise);

	// symbols not involved in the observation are supressed in "_involvedSymbolsVector"
	OneChannelOrderPerTransmitAtennaMIMOChannel::CompleteSymbolsVectorToOnlyInvolvedSymbolsVector(symbolsVector,_N,_antennasChannelOrders,_involvedSymbolsVector);

	OneChannelOrderPerTransmitAtennaMIMOChannel::WithoutZerosMatrixToWithZerosMatrix(_realEstimator->NextMatrix(observations,_involvedSymbolsVector,noiseVariance),_N,_antennasChannelOrders,_lastEstimatedChannelMatrix);

	return _lastEstimatedChannelMatrix;
}

