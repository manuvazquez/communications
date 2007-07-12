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
#ifndef MMSEDETECTOR_H
#define MMSEDETECTOR_H

#include <LinearDetector.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <Util.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>

class MMSEDetector : public LinearDetector
{
protected:
	int _nSymbolsToBeDetected;
	tMatrix _filter;

	// auxiliary variables
	tMatrix _alphabetVarianceChannelMatrixChannelMatrixTrans;
	tMatrix _Rx;
	tLongIntVector _piv;
	tVector _softEstimations;
	tRange _rNsimbolsDetected,_rAllChannelMatrixRows;

	// required for NthSymbolVariance computing
	tMatrix _channelMatrix;
	tMatrix _RxBak;
public:
    MMSEDetector(int rows, int cols, double alphabetVariance,int nSymbolsToBeDetected);

    virtual MMSEDetector * Clone();
	virtual tMatrix ComputedFilter();
    virtual tVector Detect(tVector observations, tMatrix channelMatrix, const tMatrix& noiseCovariance);
    virtual void StateStep(tVector observations) {}
	virtual double NthSymbolVariance(int n);
	virtual double NthSymbolGain(int n) const;

};

#endif
