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
#ifndef KALMANESTIMATOR_H
#define KALMANESTIMATOR_H

#include <ChannelMatrixEstimator.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <KalmanFilter.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>

class KalmanEstimator : public ChannelMatrixEstimator
{
private:
	KalmanFilter *_kalmanFilter;
	int _nChannelCoefficients;
	tMatrix _identityL;


	// variables needed for Clone
	double _ARcoefficient,_ARvariance;

	// auxiliary variables (just for efficiency's sake)
	tMatrix _F;
	tLongIntVector _piv;
	tMatrix _FtransInvNoiseCovariance,_B;
	tVector _invPredictiveCovariancePredictiveMean,_auxAuxArgExp,_auxAuxArgExpInvB,_observationsNoiseCovariance;

private:
	void FillFfromSymbolsMatrix(const tMatrix &symbolsMatrix);
public:
    KalmanEstimator(tMatrix &initialEstimation,double ARcoefficient,double ARvariance);
	~KalmanEstimator();

	tMatrix NextMatrix(const tVector &observations,const tMatrix &symbolsMatrix,double noiseVariance);
	double Likelihood(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance);
	KalmanEstimator *Clone();
};

#endif
