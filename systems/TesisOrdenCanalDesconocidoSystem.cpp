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
#include "TesisOrdenCanalDesconocidoSystem.h"

TesisOrdenCanalDesconocidoSystem::TesisOrdenCanalDesconocidoSystem()
 : ChannelOrderEstimationSystem()
{
    nSurvivors = 12;
//     nSurvivors = 1;    
    
    adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;

    forgettingFactor = 0.99;
    forgettingFactorDetector = 0.95;

    _powerProfile = new FlatPowerProfile(_L,_N,_m,1.0);

    if(adjustParticlesNumberFromSurvivors && adjustSurvivorsFromParticlesNumber)
        throw RuntimeException("adjustParticlesNumberFromSurvivors y adjustSurvivorsFromParticlesNumber no pueden ser true a la vez.");

    if(adjustParticlesNumberFromSurvivors)
    {
        nParticles = (int)pow((double)_alphabet->length(),_N*(_m-1))*nSurvivors;
        cout << "Number of particles adjusted to " << nParticles << endl;
    }

    if(adjustSurvivorsFromParticlesNumber)
    {
        cout << "number of survivors adjusted from " << nSurvivors;
        nSurvivors = int(ceil(double(nParticles)/pow(2.0,double(_N*(_m-1)))));
        cout << " to " << nSurvivors << endl;
    }

    rmmseDetector = new RMMSEDetector(_L*(c+_d+1),_N*(_m+c+_d),_alphabet->variance(),forgettingFactorDetector,_N*(_d+1));
    rlsEstimator = new RLSEstimator(_powerProfile->means(),_N,forgettingFactor);

    for(uint iChannelOrder=0;iChannelOrder<_candidateChannelOrders.size();iChannelOrder++)
    {
        RLSchannelEstimators.push_back(new RLSEstimator(_channelOrderCoefficientsMeans[iChannelOrder],_N,forgettingFactor));
        kalmanChannelEstimators.push_back(new KalmanEstimator(_channelOrderCoefficientsMeans[iChannelOrder],_channelOrderCoefficientsVariances[iChannelOrder],_N,ARcoefficients,ARvariance));
        noForgetRLSchannelEstimators.push_back(new RLSEstimator(_channelOrderCoefficientsMeans[iChannelOrder],_N,1.0));

        RMMSElinearDetectors.push_back(new RMMSEDetector(_L*_candidateChannelOrders[iChannelOrder],_N*(2*_candidateChannelOrders[iChannelOrder]-1),_alphabet->variance(),forgettingFactorDetector,_N*_candidateChannelOrders[iChannelOrder]));
    }

    ResamplingCriterion resamplingCriterion(resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
    bestParticlesResamplingAlgorithm = new BestParticlesResamplingAlgorithm(resamplingCriterion);
    multinomialResamplingAlgorithm = new MultinomialResamplingAlgorithm(resamplingCriterion);

    kalmanEstimator = new KalmanEstimator(_powerProfile->means(),_powerProfile->variances(),_N,ARcoefficients,ARvariance);

    // USIS2SIS transition criterion(s)
    USISmaximumProbabilityCriterion = new MaximumProbabilityCriterion(0.8);
//     USISuniformRelatedCriterion = new UniformRelatedCriterion(2.0);

    channelOrderEstimator= new APPbasedChannelOrderEstimator(_N,_candidateChannelOrders);
}


TesisOrdenCanalDesconocidoSystem::~TesisOrdenCanalDesconocidoSystem()
{
    delete _powerProfile;

    delete rmmseDetector;

    delete rlsEstimator;
    for(uint iChannelOrder=0;iChannelOrder<_candidateChannelOrders.size();iChannelOrder++)
    {
        delete RLSchannelEstimators[iChannelOrder];
        delete kalmanChannelEstimators[iChannelOrder];
        delete noForgetRLSchannelEstimators[iChannelOrder];
        delete RMMSElinearDetectors[iChannelOrder];
    }

    delete withoutReplacementResamplingAlgorithm;
    delete bestParticlesResamplingAlgorithm;
    delete multinomialResamplingAlgorithm;

    delete kalmanEstimator;

    delete channelOrderEstimator;
    delete USISmaximumProbabilityCriterion;
//  delete USISuniformRelatedCriterion;
}

void TesisOrdenCanalDesconocidoSystem::addAlgorithms()
{
    ChannelOrderEstimationSystem::addAlgorithms();

	// we gotta make sure the channel has memory...
	assert(_iTrueChannelOrder!=0);
	
    _algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS (subestimado)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m-1,RLSchannelEstimators[_iTrueChannelOrder-1],RMMSElinearDetectors[_iTrueChannelOrder-1],_preamble,c,_d-1,_d-1,nParticles,multinomialResamplingAlgorithm,_channelOrderCoefficientsMeans[_iTrueChannelOrder-1],_channelOrderCoefficientsVariances[_iTrueChannelOrder-1],ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    _algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,RLSchannelEstimators[_iTrueChannelOrder],RMMSElinearDetectors[_iTrueChannelOrder],_preamble,c,_d,_d,nParticles,multinomialResamplingAlgorithm,_powerProfile->means(),_powerProfile->variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

	assert((_iTrueChannelOrder+1)<_candidateChannelOrders.size());
    _algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS (sobreestimado)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m+1,RLSchannelEstimators[_iTrueChannelOrder+1],RMMSElinearDetectors[_iTrueChannelOrder+1],_preamble,c,_d+1,_d+1,nParticles,multinomialResamplingAlgorithm,_channelOrderCoefficientsMeans[_iTrueChannelOrder+1],_channelOrderCoefficientsVariances[_iTrueChannelOrder+1],ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

// // 	assert((_iTrueChannelOrder+2)<_candidateChannelOrders.size());
//     _algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS (sobreestimado por 2)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m+2,RLSchannelEstimators[_iTrueChannelOrder+2],RMMSElinearDetectors[_iTrueChannelOrder+2],_preamble,c,_d+2,_d+2,nParticles,multinomialResamplingAlgorithm,_channelOrderCoefficientsMeans[_iTrueChannelOrder+2],_channelOrderCoefficientsVariances[_iTrueChannelOrder+2],ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    _algorithms.push_back(new CMEBasedAlgorithm("CME based algorithm (RLS no forget)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,noForgetRLSchannelEstimators,_preamble,_preamble.cols(),_symbols));

    _algorithms.push_back(new TimeVaryingChannelCMEbasedAlgorithm("TimeVaryingChannelCMEbasedAlgorithm",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,RLSchannelEstimators,_preamble,_preamble.cols(),_symbols));

    // -----------------------------------------------
// //     algorithms.push_back(new LinearFilterBasedCMEapplyingAlgorithm("Linear Filter with CME based order estimation", *alphabet,L,L,N,iLastSymbolVectorToBeDetected,kalmanChannelEstimators,preamble, MMSEsmallLinearDetectors, ARcoefficients[0], true));

// //     algorithms.push_back(new KnownSymbolsCMEapplyingAlgorithm("KnownSymbolsCMEapplyingAlgorithm",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,noForgetRLSchannelEstimators,preamble,symbols));
    // -----------------------------------------------

//     algorithms.push_back(new ISIS("ISIS",*alphabet,L,L,N,iLastSymbolVectorToBeDetected,kalmanChannelEstimators,preamble,preamble.cols(),d,nParticles,multinomialResamplingAlgorithm));

    _algorithms.push_back(new USIS("USIS",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,RLSchannelEstimators,RMMSElinearDetectors,_preamble,_preamble.cols(),_d,nParticles,multinomialResamplingAlgorithm,channelOrderEstimator,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

// //     _algorithms.push_back(new USIS2SISAlgorithm("USIS2SIS",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,RLSchannelEstimators,RMMSElinearDetectors,_preamble,_preamble.cols(),_d,nParticles,multinomialResamplingAlgorithm,channelOrderEstimator,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,USISmaximumProbabilityCriterion));

    _algorithms.push_back(new MLSDmAlgorithm("MLSD-m (RLS)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,RLSchannelEstimators,_preamble,_preamble.cols(),_d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    _algorithms.push_back(new MLSDmAlgorithm("MLSD-m (Kalman)",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,kalmanChannelEstimators,_preamble,_preamble.cols(),_d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    _algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,_m,kalmanEstimator,_preamble,_d,_iLastSymbolVectorToBeDetected+_d,nSurvivors));

    _algorithms.push_back(new ViterbiAlgorithm("Viterbi",*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,*(dynamic_cast<StillMemoryMIMOChannel *> (_channel)),_preamble,_d));
}

void TesisOrdenCanalDesconocidoSystem::beforeEndingFrame()
{
    ChannelOrderEstimationSystem::beforeEndingFrame();
    Util::scalarToOctaveFileStream(nSurvivors,"nSurvivors",_f);
    Util::scalarToOctaveFileStream(forgettingFactor,"forgettingFactor",_f);
    Util::scalarToOctaveFileStream(forgettingFactorDetector,"forgettingFactorDetector",_f);
}

void TesisOrdenCanalDesconocidoSystem::buildSystemSpecificVariables()
{
	_noise = new PowerProfileDependentNoise(_alphabet->variance(),_L,_channel->length(),*_powerProfile);
}

void TesisOrdenCanalDesconocidoSystem::saveFrameResults()
{
    ChannelOrderEstimationSystem::saveFrameResults();
    Octave::toOctaveFileStream(nSurvivors,"nSurvivors",_f);
    Octave::toOctaveFileStream(forgettingFactor,"forgettingFactor",_f);
    Octave::toOctaveFileStream(forgettingFactorDetector,"forgettingFactorDetector",_f);
}
