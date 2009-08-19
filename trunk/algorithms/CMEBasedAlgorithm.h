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
#ifndef CMEBASEDALGORITHM_H
#define CMEBASEDALGORITHM_H

#include <UnknownChannelOrderAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <RLSEstimator.h>
#include <Eigen/LU> 

class CMEBasedAlgorithm : public UnknownChannelOrderAlgorithm
{
protected:
    MatrixXd _symbolVectors;
public:
    CMEBasedAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, const tMatrix &symbolVectors);

    virtual void run(tMatrix observations,vector<double> noiseVariances)
    {
        run(Util::lapack2eigen(observations),noiseVariances);
    }
    virtual void run(MatrixXd observations,vector<double> noiseVariances);
    
    virtual void run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
    {
        run(Util::lapack2eigen(observations),noiseVariances,Util::lapack2eigen(trainingSequence));
    }
    virtual void run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence);

    virtual tMatrix getDetectedSymbolVectors()
    {
        return Util::eigen2lapack(getDetectedSymbolVectors_eigen());
    }
    virtual MatrixXd getDetectedSymbolVectors_eigen();
    
    virtual vector<tMatrix> getEstimatedChannelMatrices()
    {
        return Util::eigen2lapack(getEstimatedChannelMatrices_eigen());
    }
    virtual vector<MatrixXd> getEstimatedChannelMatrices_eigen();      
};

#endif
