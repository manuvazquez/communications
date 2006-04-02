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
#ifndef KALMANFILTER_H
#define KALMANFILTER_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <exceptions.h>
#include <Util.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>

class KalmanFilter{
private:
	tMatrix _R, _Rtranspose, _stateEquationCovariance;
	int _nElementsToEstimate;
	tVector _predictiveMean,_filteredMean;
	tMatrix _predictiveCovariance,_filteredCovariance;
	int _observationVectorLength;

	// auxiliar variables
	tMatrix _predictiveCovarianceFtrans,_auxMatrix,_KalmanGain,_FpredictiveCovariance;
	tMatrix _KalmanGainFpredictiveCovariance,_predictiveCovarianceAux;
	tMatrix _RfilteredCovariance,_RfilteredCovarianceRtrans;
	tVector _auxVector,_KalmanGainByNotPredicted;
	tLongIntVector _piv;
public:
    KalmanFilter(tMatrix R,tMatrix stateEquationCovariance,tVector initialMean, tMatrix initialCovariance,int observationVectorLength);

	void Step(tMatrix F,tVector observation, tMatrix observationEquationCovariance);
	tVector PredictiveMean() { return _predictiveMean;}
	tVector FilteredMean() {return _filteredMean;}
	tMatrix PredictiveCovariance() {return _predictiveCovariance;}
	tMatrix FilteredCovariance() {return _filteredCovariance;}
};

#endif
