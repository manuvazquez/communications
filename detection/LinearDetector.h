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
#ifndef LINEARDETECTOR_H
#define LINEARDETECTOR_H

/**
    @author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <Util.h>

class LinearDetector{
protected:
    const uint _channelMatrixRows, _channelMatrixCols;
    const double _alphabetVariance;
	MatrixXd _filter;
public:
    LinearDetector(uint rows,uint cols,double alphabetVariance);
    virtual ~LinearDetector() {}
    
    virtual void stateStep(VectorXd observations) = 0;
    virtual VectorXd detect(VectorXd observations,MatrixXd channelMatrix,const MatrixXd &noiseCovariance) = 0;
	virtual MatrixXd computedFilter() const { return _filter;}

    /**
     *    Computes the variance related to the soft estimation provided for the n-th symbol. It must NEVER be called before a call to "detect"
     * @param n
     * @return
     */
    virtual double nthSymbolVariance(uint n,double noiseVariance) = 0;

    virtual double nthSymbolGain(uint n) const { return 1.0;}
    uint channelMatrixcols() { return _channelMatrixCols;}
    virtual LinearDetector *clone() = 0;
    
    void stateStepsFromObservationsSequence(const MatrixXd &observations,uint smoothingLag,uint iFrom,uint iTo);
};

#endif
