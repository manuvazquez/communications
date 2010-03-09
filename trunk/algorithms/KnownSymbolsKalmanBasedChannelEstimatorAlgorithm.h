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
#ifndef KNOWNSYMBOLSKALMANBASEDCHANNELESTIMATORALGORITHM_H
#define KNOWNSYMBOLSKALMANBASEDCHANNELESTIMATORALGORITHM_H

#include <KnownChannelOrderAlgorithm.h>

/**
@author Manu
*/

#include <KalmanEstimator.h>

class KnownSymbolsKalmanBasedChannelEstimatorAlgorithm : public KnownChannelOrderAlgorithm
{
protected:
    MatrixXd _symbolVectors;
    vector<MatrixXd> _estimatedChannelMatrices;
public:

    /**
     * It must receive a KalmanEstimator (or a wrapped one)
     * @param name
     * @param alphabet
     * @param channelEstimator
     * @param preamble
     * @param symbolVectors includes the preamble
     * @return
     */
    KnownSymbolsKalmanBasedChannelEstimatorAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected,int m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble,const MatrixXd &symbolVectors);

    virtual void run(MatrixXd observations,vector<double> noiseVariances);
    virtual void run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence);

    virtual MatrixXd getDetectedSymbolVectors();
    
    virtual vector<MatrixXd> getEstimatedChannelMatrices();    
};

#endif
