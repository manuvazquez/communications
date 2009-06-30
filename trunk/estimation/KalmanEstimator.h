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
#include <StatUtil.h>
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

private:
	tMatrix BuildFfromSymbolsMatrix(const tVector &symbolsVector);
public:
    KalmanEstimator(const tMatrix &initialEstimation,int N,double ARcoefficient,double ARvariance);
    KalmanEstimator(const tMatrix &initialEstimation,const tMatrix &variances,int N,double ARcoefficient,double ARvariance);
	KalmanEstimator(const KalmanEstimator &kalmanEstimator);
	~KalmanEstimator();

	virtual tMatrix nextMatrix(const tVector &observations,const tMatrix &symbolsMatrix,double noiseVariance);
	double likelihood(const tVector &observations,const tMatrix symbolsMatrix,double noiseVariance);
	virtual KalmanEstimator *Clone() const;
	tMatrix SampleFromPredictive();
    virtual void setFirstEstimatedChannelMatrix(const tMatrix &matrix);
};

#endif
