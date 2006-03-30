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
#include "KalmanFilter.h"

KalmanFilter::KalmanFilter(tMatrix R,tMatrix stateEquationCovariance,tVector initialMean, tMatrix initialCovariance):
_nElementsToEstimate(initialMean.size()),_R(R),_stateEquationCovariance(stateEquationCovariance),_filteredMean(initialMean),_filteredCovariance(initialCovariance)
{
	if(R.rows()!=_nElementsToEstimate || _nElementsToEstimate!=R.cols())
		throw RuntimeException("Matrices R and F dimensions are not coherent with those of the vector to be estimated.");
	
}


KalmanFilter::~KalmanFilter()
{
}


