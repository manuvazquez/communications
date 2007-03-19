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
#ifndef RMMSEDETECTOR_H
#define RMMSEDETECTOR_H

#include <LinearDetector.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <exceptions.h>
#include <Util.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>

class RMMSEDetector : public LinearDetector
{
protected:
	double _forgettingFactor,_invForgettingFactor;
	int _nSymbolsToBeDetected;
	double _alphaPowerSumNow,_alphaPowerSumPrevious;
	double _alphaPower,_alphaPowerSumFactor;

	tVector _g;
	tMatrix _invRtilde;
	tMatrix _filter;

	// auxiliary
// 	tMatrix _identityL,_gObservations;
// 	tMatrix _identityMinusgObservations;
    tMatrix _auxInvRtilde;
	tMatrix _E,_varianceInvRtildeChannelMatrix;
public:
    RMMSEDetector(int rows, int cols,double alphabetVariance,double forgettingFactor,int nSymbolsToBeDetected);

    tVector Detect(tVector observations, tMatrix channelMatrix,const tMatrix &noiseCovariance);
	RMMSEDetector *Clone();
	void StateStep(tVector observations);

	tMatrix ComputedFilter() { return _filter;}

};

#endif
