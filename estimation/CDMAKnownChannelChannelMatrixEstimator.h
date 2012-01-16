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
#ifndef CDMAKNOWNCHANNELCHANNELMATRIXESTIMATOR_H
#define CDMAKNOWNCHANNELCHANNELMATRIXESTIMATOR_H

#include <KnownChannelChannelMatrixEstimator.h>

/**
	@author Manu <manu@rustneversleeps>
*/
class CDMAKnownChannelChannelMatrixEstimator : public KnownChannelChannelMatrixEstimator
{
protected:
    MatrixXd _spreadingCodes;
    
    virtual double likelihood(const VectorXd &observations,const MatrixXd symbolsMatrix,double noiseVariance); // eigen
public:
    CDMAKnownChannelChannelMatrixEstimator(const MIMOChannel *channel, uint iFirstChannelMatrix, uint N, const MatrixXd &spreadingCodes);

    virtual CDMAKnownChannelChannelMatrixEstimator *clone() const;
	
	virtual MatrixXd nextMatrix(const VectorXd& observations, const MatrixXd& symbolsMatrix, double noiseVariance) { return _spreadingCodes * KnownChannelChannelMatrixEstimator::nextMatrix(observations,symbolsMatrix,noiseVariance).asDiagonal();}
	
	virtual MatrixXd predictedMatrix() const { return _spreadingCodes * KnownChannelChannelMatrixEstimator::predictedMatrix().asDiagonal();}
	
	virtual MatrixXd lastEstimatedChannelMatrix() const { return _spreadingCodes * KnownChannelChannelMatrixEstimator::lastEstimatedChannelCoefficientsMatrix().asDiagonal();}
};

#endif
