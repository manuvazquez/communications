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
#include "MultipleChannelEstimatorsPerParticleSMCAlgorithm.h"

MultipleChannelEstimatorsPerParticleSMCAlgorithm::MultipleChannelEstimatorsPerParticleSMCAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation,int smoothingLag,int nParticles,ResamplingAlgorithm *resamplingAlgorithm): UnknownChannelOrderAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation)
//variables initialization
,_resamplingAlgorithm(resamplingAlgorithm),_d(smoothingLag),_allSymbolsRows(0,_N-1)
{
    // at first, we assume that all symbol vectors from the preamble need to be processed
    _startDetectionTime = _preamble.cols();
}


MultipleChannelEstimatorsPerParticleSMCAlgorithm::~MultipleChannelEstimatorsPerParticleSMCAlgorithm()
{
}

void MultipleChannelEstimatorsPerParticleSMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    int nObservations = observations.cols();

    if(nObservations<_startDetectionTime+_maxOrder)
        throw RuntimeException("MultipleChannelEstimatorsPerParticleSMCAlgorithm::Run: Not enough observations.");

    this->InitializeParticles();

    this->Process(observations,noiseVariances);
}

void MultipleChannelEstimatorsPerParticleSMCAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
    if(observations.rows()!=_L || trainingSequence.rows()!=_N)
        throw RuntimeException("MultipleChannelEstimatorsPerParticleSMCAlgorithm::Run: Observations matrix or training sequence dimensions are wrong.");

    int iParticle;
    uint j;
    uint iChannelOrder;

    // to process the training sequence, we need both the preamble and the symbol vectors related to it
    tMatrix preambleTrainingSequence = Util::Append(_preamble,trainingSequence);


    tRange rSymbolVectorsTrainingSequece(0,preambleTrainingSequence.cols()-1);

    vector<vector<tMatrix> > trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

    this->InitializeParticles();

    for(iParticle=0;iParticle<GetParticleFilterPointer()->Nparticles();iParticle++)
    {
        ParticleWithChannelEstimation *processedParticle = dynamic_cast <ParticleWithChannelEstimation *> (GetParticleFilterPointer()->GetParticle(iParticle));

        for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        {
            //the channel estimation given by the training sequence is copied into each particle...
            for(j=0;j<trainingSequenceChannelMatrices[iChannelOrder].size();j++)
            {
                processedParticle->SetChannelMatrix(iChannelOrder,_preamble.cols()+j,trainingSequenceChannelMatrices[iChannelOrder][j]);
            }
        }

        //... the symbols are considered detected...
        processedParticle->SetSymbolVectors(rSymbolVectorsTrainingSequece,preambleTrainingSequence);
    }

    // the Process method must start in
    _startDetectionTime += trainingSequence.cols();

    this->Process(observations,noiseVariances);
}

tMatrix MultipleChannelEstimatorsPerParticleSMCAlgorithm::GetDetectedSymbolVectors()
{
//     int iBestParticle;
//     Util::Max(GetParticleFilterPointer()->GetWeightsVector(),iBestParticle);


    return (GetParticleFilterPointer()->GetBestParticle()->GetAllSymbolVectors())(_allSymbolsRows,tRange(_preamble.cols(),_K-1));
}

vector<tMatrix> MultipleChannelEstimatorsPerParticleSMCAlgorithm::GetEstimatedChannelMatrices()
{
    vector<tMatrix> channelMatrices;
    channelMatrices.reserve(_K-_preamble.cols());

    // best particle is chosen
	int iBestParticle = GetParticleFilterPointer()->iBestParticle();
//     Util::Max(GetParticleFilterPointer()->GetWeightsVector(),iBestParticle);

	int iBestChannelOrder = BestChannelOrderIndex(iBestParticle);

    for(int i=_preamble.cols();i<_K;i++)
        channelMatrices.push_back((GetParticleFilterPointer()->GetParticle(iBestParticle))->GetChannelMatrix(iBestChannelOrder,i));

    return channelMatrices;
}
