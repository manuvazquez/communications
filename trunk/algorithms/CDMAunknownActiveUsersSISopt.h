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
#ifndef CDMAUNKNOWNACTIVEUSERSSISOPT_H
#define CDMAUNKNOWNACTIVEUSERSSISOPT_H

#include <SMCAlgorithm.h>

/**
It implements an (optimal) algorithm that aims to detect the active users in a SISO CDMA system along with the transmitted data

    @author Manu <manu@rustneversleeps>
*/
class CDMAunknownActiveUsersSISopt : public SMCAlgorithm
{
public:
    CDMAunknownActiveUsersSISopt(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm, const tMatrix& channelMatrixMean, const tMatrix& channelMatrixVariances);

    ~CDMAunknownActiveUsersSISopt();

    tMatrix getDetectedSymbolVectors();
    vector< tMatrix > GetEstimatedChannelMatrices();

protected:
    virtual void BeforeInitializingParticles(const tMatrix& observations, const tMatrix& trainingSequence);
    virtual void InitializeParticles();
    virtual void Process(const tMatrix& observations, vector< double > noiseVariances);

};

#endif
