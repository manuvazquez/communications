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
#ifndef UTRELLISSEARCHALGORITHM_H
#define UTRELLISSEARCHALGORITHM_H

#include <ChannelOrderEstimatorSMCAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/

#include <ParticleWithChannelEstimationAndChannelOrderAPP.h>

class UTrellisSearchAlgorithm : public ChannelOrderEstimatorSMCAlgorithm
{
protected:
	ParticleFilter *_particleFilter;
	double _ARcoefficient,_samplingVariance,_ARprocessVariance;
    vector<int> _particlesBestChannelOrders;

	vector<vector<tMatrix> > ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence);
    int BestChannelOrderIndex(int iBestParticle);
public:
    UTrellisSearchAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm,double ARcoefficient,double samplingVariance,double ARprocessVariance);

    ~UTrellisSearchAlgorithm();

    virtual ParticleFilter* GetParticleFilterPointer() {return _particleFilter;}
    virtual void InitializeParticles();
    virtual void Process(const tMatrix& observations, vector< double > noiseVariances);

};

#endif
