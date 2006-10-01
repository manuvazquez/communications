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
#include "LinearFilterBasedISIRAlgorithm.h"

LinearFilterBasedISIRAlgorithm::LinearFilterBasedISIRAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators,vector<LinearDetector *> linearDetectors, tMatrix preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm): MultipleChannelEstimatorsPerParticleSMCAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm),_allObservationRows(0,_L-1),_linearDetectors(linearDetectors.size())
{
    if(linearDetectors.size()!=_candidateOrders.size())
        throw RuntimeException("LinearFilterBasedISIRAlgorithm::LinearFilterBasedISIRAlgorithm: nº of detectors and number of channel matrix estimators (and candidate orders) are different.");

    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _linearDetectors[iChannelOrder] = linearDetectors[iChannelOrder]->Clone();
}


LinearFilterBasedISIRAlgorithm::~LinearFilterBasedISIRAlgorithm()
{
    for(int iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        delete _linearDetectors[iChannelOrder];
}

vector<vector<tMatrix> > LinearFilterBasedISIRAlgorithm::ProcessTrainingSequence(const tMatrix &observations,vector<double> noiseVariances,tMatrix trainingSequence)
{
    int iChannelOrder;

    for(int i=_iFirstObservation;i<_iFirstObservation+trainingSequence.cols();i++)
    {
        for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        {
            // the observations from i to i+d are stacked
            tRange smoothingRange(i,i+_candidateOrders[iChannelOrder]-1);
            tVector stackedObservationsVector = Util::ToVector(observations(_allObservationRows,smoothingRange),columnwise);

            _linearDetectors[iChannelOrder]->StateStep(stackedObservationsVector);
        }
    }

    return UnknownChannelOrderAlgorithm::ProcessTrainingSequence(observations,noiseVariances,trainingSequence);
}

ParticleFilter* LinearFilterBasedISIRAlgorithm::GetParticleFilterPointer()
{
}

void LinearFilterBasedISIRAlgorithm::InitializeParticles()
{
}

void LinearFilterBasedISIRAlgorithm::Process(const tMatrix& observations, vector< double > noiseVariances)
{
}

