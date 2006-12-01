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
#ifndef LMSESTIMATOR_H
#define LMSESTIMATOR_H

#include <ChannelMatrixEstimator.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <Util.h>
#include <lapackpp/gmd.h>
// #include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
// #include <lapackpp/blas3pp.h>
// #include <lapackpp/laslv.h>
// #include <lapackpp/lavli.h>

class LMSEstimator : public ChannelMatrixEstimator
{
protected:
	double _mu;
	tVector _predictedObservations,_error;
	tMatrix _deltaMatrix;
public:
    LMSEstimator(const tMatrix &initialEstimation,double mu);

    virtual LMSEstimator* Clone();
    virtual tMatrix NextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance);

};

#endif
