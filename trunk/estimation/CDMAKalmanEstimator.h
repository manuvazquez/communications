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
#ifndef CDMAKALMANESTIMATOR_H
#define CDMAKALMANESTIMATOR_H

#include <ChannelMatrixEstimator.h>
#include <KalmanFilter.h>
#include <Util.h>

/**
It implements a channel matrix estimator for a Multiuser CDMA autoregressive channel

	@author Manu <manu@rustneversleeps>
*/
class CDMAKalmanEstimator : public ChannelMatrixEstimator
{
private:
    tMatrix BuildFfromSymbolsMatrix(const tVector &symbolsVector);
protected:
    uint _stateVectorLength,_nOutputs;
    KalmanFilter *_kalmanFilter;
    tMatrix _spreadingCodes;
public:
    CDMAKalmanEstimator(tMatrix initialEstimation, const tMatrix &spreadingCodes, const vector<double> &ARprocCoeffs, double ARprocVar);
    
    CDMAKalmanEstimator(int N);

    ~CDMAKalmanEstimator();

    virtual ChannelMatrixEstimator* Clone() const;
    virtual double likelihood(const tVector& observations, const tMatrix symbolsMatrix, double noiseVariance);
    virtual tMatrix nextMatrix(const tVector& observations, const tMatrix& symbolsMatrix, double noiseVariance);

};

#endif
