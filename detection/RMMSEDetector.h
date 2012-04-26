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

class RMMSEDetector : public LinearDetector
{
protected:
    double _forgettingFactor,_invForgettingFactor;
    uint _nSymbolsToBeDetected;
    double _alphaPowerSumNow,_alphaPowerSumPrevious;
    double _alphaPower,_alphaPowerSumFactor;

    VectorXd _g;
    MatrixXd _invRtilde,_filter;

    MatrixXd _E;

    // required for nthSymbolVariance computing
    MatrixXd _channelMatrix,_alphabetVarianceChannelMatrixChannelMatrixTransPlusNoiseCovariance;
public:
    RMMSEDetector(uint rows, uint cols,double alphabetVariance,double forgettingFactor,uint nSymbolsToBeDetected);

    virtual VectorXd detect(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance);

    RMMSEDetector *clone();
    void stateStep(VectorXd observations);

    double nthSymbolVariance(uint n,double noiseVariance);
    MatrixXd computedFilter() const { return _filter;}

};

#endif
