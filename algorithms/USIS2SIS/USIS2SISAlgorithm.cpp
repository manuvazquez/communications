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
#include "USIS2SISAlgorithm.h"

// #define DEBUG

USIS2SISAlgorithm::USIS2SISAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, vector< LinearDetector * > linearDetectors, MatrixXd preamble, uint iFirstObservation, uint smoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm, ChannelOrderEstimator* channelOrderEstimator, double ARcoefficient, double samplingVariance, double ARprocessVariance, TransitionCriterion *transitionCriterion): USIS(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, linearDetectors, preamble, iFirstObservation, smoothingLag, nParticles, resamplingAlgorithm, channelOrderEstimator, ARcoefficient, samplingVariance, ARprocessVariance),_transitionCriterion(transitionCriterion)
{
}

void USIS2SISAlgorithm::beforeResamplingProcess(uint iProcessedObservation, const MatrixXd& observations, const vector<double> &noiseVariances)
{
    VectorXd _weightedChannelOrderAPPs = VectorXd::Zero(_candidateOrders.size());

    for(uint iParticle=0;iParticle<_particleFilter.capacity();iParticle++)
    {
        ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *processedParticle = dynamic_cast <ParticleWithChannelEstimationAndLinearDetectionAndChannelOrderEstimation *>(_particleFilter.getParticle(iParticle));

        VectorXd particleChannelOrderAPPs = processedParticle->getChannelOrderEstimator()->getChannelOrderAPPsVector();
        _weightedChannelOrderAPPs += processedParticle->getWeight()*particleChannelOrderAPPs;
    }

    // the maximum probability is obtained
    uint iMax;
    _weightedChannelOrderAPPs.maxCoeff(&iMax);

    // if the transition criterion is satisfied
    if(_transitionCriterion->makeTransition(_weightedChannelOrderAPPs))
    {
        LinearFilterBasedSMCAlgorithm knownChannelOrderAlgorithm(_name,_alphabet,_nOutputs,_nOutputs,_nInputs,_iLastSymbolVectorToBeDetected,_candidateOrders[iMax],_preamble,_candidateOrders[iMax]-1,&_particleFilter,_resamplingAlgorithm,_ARcoefficient,_samplingVariance,_ARprocessVariance);

        knownChannelOrderAlgorithm.setEstimatorIndex(iMax);
        knownChannelOrderAlgorithm.runFrom(iProcessedObservation,observations,noiseVariances);
        _processDoneExternally = true;

#ifdef DEBUG
		cout << "USIS2SIS: transition criterion holds (iProcessedObservation = " << iProcessedObservation << " )" << endl;
#endif

        // the APP of the selected channel order is set to 1.0
        _channelOrderAPPs[0].row(iMax).segment(iProcessedObservation,_iLastSymbolVectorToBeDetected-iProcessedObservation).setConstant(1.0);

        return;
    }
}
