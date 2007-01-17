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
#include "OneChannelOrderPerTransmitAtennaRLSEstimator.h"

// #define DEBUG

OneChannelOrderPerTransmitAtennaRLSEstimator::OneChannelOrderPerTransmitAtennaRLSEstimator(const tMatrix& initialEstimation,int N,double forgettingFactor,const vector<int> &antennasChannelOrders): RLSEstimator(N,forgettingFactor),_antennasChannelOrders(antennasChannelOrders)
{
	if(antennasChannelOrders.size()!=_N)
		throw RuntimeException("OneChannelOrderPerTransmitAtennaRLSEstimator::OneChannelOrderPerTransmitAtennaRLSEstimator: antennasChannelOrders size is wrong.");


	// initialization of variables declared in ChannelMatrixEstimator
	_L = initialEstimation.rows();
	_Nm = Util::Sum(_antennasChannelOrders);
	if((initialEstimation.cols() % _N)!=0)
		throw RuntimeException("OneChannelOrderPerTransmitAtennaRLSEstimator::OneChannelOrderPerTransmitAtennaRLSEstimator: the length of the antennas channel orders vector is not coherent with the initial estimation channel matrix dimensions");

	_lastEstimatedChannelMatrix = tMatrix(_L,_Nm);

	// initialization of variables declared in RLSEstimator
    _invRtilde = LaGenMatDouble::eye(_Nm);
    _pTilde = LaGenMatDouble::zeros(_L,_Nm);
    _invForgettingFactorSymbolsVectorInvRtilde = tVector(_Nm);
    _g = tVector(_Nm);
    _invForgettingFactorInvRtildeSymbolsVector = tVector(_Nm);
    _invForgettingFactorInvRtildeSymbolsVectorg = tMatrix(_Nm,_Nm);
    _observationsSymbolsVector = tMatrix(_L,_Nm);

	// initialization of the wrapper variables
    _involvedSymbolsVector = tVector(_Nm);
//     _withZerosMatrix = tMatrix(initialEstimation.rows(),initialEstimation.cols());
    _withZerosMatrix = initialEstimation;
    _withoutZerosMatrix = tMatrix(_L,_Nm);

    #ifdef DEBUG
    	cout << "Antes de convertir initialEstimation" << endl;
    #endif

	// the initial estimation received is cut according to the orders of the transmit antennas
	OneChannelOrderPerTransmitAtennaMIMOChannel::WithZerosMatrixToWithoutZerosMatrix(initialEstimation,_N,_antennasChannelOrders,_lastEstimatedChannelMatrix);

    #ifdef DEBUG
    	cout << "Despues de convertir initialEstimation" << endl;
    	cout << "_lastEstimatedChannelMatrix es" << endl << _lastEstimatedChannelMatrix << endl;
    #endif
}

ChannelMatrixEstimator *OneChannelOrderPerTransmitAtennaRLSEstimator::Clone()
{
	return new OneChannelOrderPerTransmitAtennaRLSEstimator(*this);
}

tMatrix OneChannelOrderPerTransmitAtennaRLSEstimator::LastEstimatedChannelMatrix()
{
	return _withZerosMatrix;
}

tMatrix OneChannelOrderPerTransmitAtennaRLSEstimator::NextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance)
{
	#ifdef DEBUG
		cout << "en el netxMatrix del wrapper" << endl;
	#endif
	tVector symbolsVector = Util::ToVector(symbolsMatrix,columnwise);

	// a new vector which includes only the symbols involving the observations is constructed
	OneChannelOrderPerTransmitAtennaMIMOChannel::CompleteSymbolsVectorToOnlyInvolvedSymbolsVector(symbolsVector,_N,_antennasChannelOrders,_involvedSymbolsVector);

	#ifdef DEBUG
		cout << "antes de llamar a NextMatrix" << endl;
	#endif

	// the underlying RLS estimator is called
	_withoutZerosMatrix = RLSEstimator::NextMatrixProcessing(observations, _involvedSymbolsVector, noiseVariance);

	// the returned matrix is converted
	OneChannelOrderPerTransmitAtennaMIMOChannel::WithoutZerosMatrixToWithZerosMatrix(_withoutZerosMatrix,_N,_antennasChannelOrders,_withZerosMatrix);


	#ifdef DEBUG
		cout << "antes del retorno del netxMatrix del wrapper" << endl;
		cout << "Va a devolver" << endl << _withZerosMatrix << endl;
		cout << "Filas: " << _withZerosMatrix.rows() << " Columnas: " << _withZerosMatrix.cols() << endl;
		getchar();
	#endif
    return _withZerosMatrix;
}

