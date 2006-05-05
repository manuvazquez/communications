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
	double _invForgettingFactor;
	int _nSymbolsToBeDetected;

	tVector _g;
	tMatrix _invRtilde;
	tMatrix _filter;

	// auxiliary
	tMatrix _identityL,_gObservations;
	tMatrix _identityMinusgObservations,_auxInvRtilde;
	tMatrix _E,_varianceInvRtildeChannelMatrix;
public:
    RMMSEDetector(int rows, int cols,double alphabetVariance,double forgettingFactor,int nSymbolsToBeDetected);

    tVector Detect(tVector observations, tMatrix channelMatrix);
	RMMSEDetector *Clone();

	void StateStep(tVector observations)
	{
		if(observations.size()!= _channelMatrixRows)
			throw RuntimeException("observations vector dimensions are wrong.");
	
		// _g = _invForgettingFactor*_invRtilde*observations
		Blas_Mat_Vec_Mult(_invRtilde,observations,_g,_invForgettingFactor);
	
		// _g = _g / (1 + _invForgettingFactor*observations'*_invRtilde*observations
		_g *= 1.0/(1.0 + Blas_Dot_Prod(observations,_g));
	
		// _gObservations = _g*observations
		Util::Mult(_g,observations,_gObservations);
	
		// _identityMinusgObservations = _identityL - _gObservations
		Util::Add(_identityL,_gObservations,_identityMinusgObservations,1.0,-1.0);
	
		// _auxInvRtilde = _invForgettingFactor*_identityMinusgObservations*_invRtilde
		Blas_Mat_Mat_Mult(_identityMinusgObservations,_invRtilde,_auxInvRtilde,_invForgettingFactor);

		_invRtilde = _auxInvRtilde;
	}

	tMatrix ComputedFilter() { return _filter;}

};

#endif
