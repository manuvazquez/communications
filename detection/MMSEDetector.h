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
#include <Eigen/LU>

class MMSEDetector : public LinearDetector
{
protected:
	int _nSymbolsToBeDetected,_detectionStart;
    MatrixXd _filter_eigen;   

	// required for nthSymbolVariance computing
    MatrixXd _channelMatrix_eigen,_Rx_eigen;   
public:
    MMSEDetector(int rows, int cols, double alphabetVariance,int nSymbolsToBeDetected);
    MMSEDetector(int rows, int cols, double alphabetVariance,int nSymbolsToBeDetected,int startingFrom);

    virtual MMSEDetector * clone();
    virtual MatrixXd computedFilter_eigen();
    virtual VectorXd detect(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance); // eigen
    virtual void stateStep(VectorXd observations) {}
	virtual double nthSymbolVariance(int n,double noiseVariance);
	virtual double nthSymbolGain(int n) const;

};

#endif
