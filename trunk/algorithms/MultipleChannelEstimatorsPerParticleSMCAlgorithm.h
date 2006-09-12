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
#ifndef MULTIPLECHANNELESTIMATORSPERPARTICLESMCALGORITHM_H
#define MULTIPLECHANNELESTIMATORSPERPARTICLESMCALGORITHM_H

#include <UnknownChannelOrderAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <ParticleFilterWithChannelOrder.h>
#include <ResamplingAlgorithm.h>

class MultipleChannelEstimatorsPerParticleSMCAlgorithm : public UnknownChannelOrderAlgorithm
{
protected:
    ParticleFilterWithChannelOrder _particleFilter;
    ResamplingAlgorithm *_resamplingAlgorithm;
    int _d,_startDetectionObservation,_startDetectionSymbolVector;
	int _startDetectionTime;
    tRange _allSymbolsRows;
    vector<int> _nParticlesPerChannelOrder;
//     tMatrix _observations;

    virtual void InitializeParticles();
    virtual void Process(const tMatrix &observations,vector<double> noiseVariances) = 0;
    vector<vector<int> > GetIndexesOfChannelOrders();
    int BestParticle();
public:
    MultipleChannelEstimatorsPerParticleSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation,int smoothingLag,int nParticles,ResamplingAlgorithm *resamplingAlgorithm);

    ~MultipleChannelEstimatorsPerParticleSMCAlgorithm();

    void Run(tMatrix observations,vector<double> noiseVariances);
    void Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);
    tMatrix GetDetectedSymbolVectors();
    vector<tMatrix> GetEstimatedChannelMatrices();

};

#endif