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
    adjustSurvivorsFromParticlesNumber = false;
    adjustParticlesNumberFromSurvivors = true;

    forgettingFactor = 0.99;
    forgettingFactorDetector = 0.95;

    powerProfile = new FlatPowerProfile(L,N,m,1.0);

    if(adjustParticlesNumberFromSurvivors && adjustSurvivorsFromParticlesNumber)
        throw RuntimeException("adjustParticlesNumberFromSurvivors y adjustSurvivorsFromParticlesNumber no pueden ser true a la vez.");

    if(adjustParticlesNumberFromSurvivors)
    {
        nParticles = (int)pow((double)alphabet->Length(),N*(m-1))*nSurvivors;
        cout << "Number of particles adjusted to " << nParticles << endl;
    }

    if(adjustSurvivorsFromParticlesNumber)
    {
        cout << "number of survivors adjusted from " << nSurvivors;
        nSurvivors = int(ceil(double(nParticles)/pow(2.0,double(N*(m-1)))));
        cout << " to " << nSurvivors << endl;
    }

    rmmseDetector = new RMMSEDetector(L*(c+d+1),N*(m+c+d),alphabet->Variance(),forgettingFactorDetector,N*(d+1));
    rlsEstimator = new RLSEstimator(powerProfile->Means(),N,forgettingFactor);

    for(uint iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
    {
        RLSchannelEstimators.push_back(new RLSEstimator(channelOrderCoefficientsMeans[iChannelOrder],N,forgettingFactor));
        kalmanChannelEstimators.push_back(new KalmanEstimator(channelOrderCoefficientsMeans[iChannelOrder],channelOrderCoefficientsVariances[iChannelOrder],N,ARcoefficients[0],ARvariance));
        noForgetRLSchannelEstimators.push_back(new RLSEstimator(channelOrderCoefficientsMeans[iChannelOrder],N,1.0));

        RMMSElinearDetectors.push_back(new RMMSEDetector(L*candidateChannelOrders[iChannelOrder],N*(2*candidateChannelOrders[iChannelOrder]-1),alphabet->Variance(),forgettingFactorDetector,N*candidateChannelOrders[iChannelOrder]));
    }

    ResamplingCriterion resamplingCriterion(resamplingRatio);
    withoutReplacementResamplingAlgorithm = new WithoutReplacementResamplingAlgorithm(resamplingCriterion);
    bestParticlesResamplingAlgorithm = new BestParticlesResamplingAlgorithm(resamplingCriterion);
    multinomialResamplingAlgorithm = new MultinomialResamplingAlgorithm(resamplingCriterion);

    kalmanEstimator = new KalmanEstimator(powerProfile->Means(),powerProfile->Variances(),N,ARcoefficients[0],ARvariance);

    // USIS2SIS transition criterion(s)
    USISmaximumProbabilityCriterion = new MaximumProbabilityCriterion(0.8);
//     USISuniformRelatedCriterion = new UniformRelatedCriterion(2.0);

    channelOrderEstimator= new APPbasedChannelOrderEstimator(N,candidateChannelOrders);
}


TesisOrdenCanalDesconocidoSystem::~TesisOrdenCanalDesconocidoSystem()
{
    delete powerProfile;

    delete rmmseDetector;

    delete rlsEstimator;
    for(uint iChannelOrder=0;iChannelOrder<candidateChannelOrders.size();iChannelOrder++)
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

void TesisOrdenCanalDesconocidoSystem::AddAlgorithms()
{
    ChannelOrderEstimationSystem::AddAlgorithms();

    algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS (subestimado)",*alphabet,L,N,lastSymbolVectorInstant,m-1,RLSchannelEstimators[iTrueChannelOrder-1],RMMSElinearDetectors[iTrueChannelOrder-1],preamble,c,d-1,d-1,nParticles,multinomialResamplingAlgorithm,channelOrderCoefficientsMeans[iTrueChannelOrder-1],channelOrderCoefficientsVariances[iTrueChannelOrder-1],ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS",*alphabet,L,N,lastSymbolVectorInstant,m,RLSchannelEstimators[iTrueChannelOrder],RMMSElinearDetectors[iTrueChannelOrder],preamble,c,d,d,nParticles,multinomialResamplingAlgorithm,powerProfile->Means(),powerProfile->Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS (sobreestimado)",*alphabet,L,N,lastSymbolVectorInstant,m+1,RLSchannelEstimators[iTrueChannelOrder+1],RMMSElinearDetectors[iTrueChannelOrder+1],preamble,c,d+1,d+1,nParticles,multinomialResamplingAlgorithm,channelOrderCoefficientsMeans[iTrueChannelOrder+1],channelOrderCoefficientsVariances[iTrueChannelOrder+1],ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    algorithms.push_back(new LinearFilterBasedSMCAlgorithm("RLS-D-SIS (sobreestimado por 2)",*alphabet,L,N,lastSymbolVectorInstant,m+2,RLSchannelEstimators[iTrueChannelOrder+2],RMMSElinearDetectors[iTrueChannelOrder+2],preamble,c,d+2,d+2,nParticles,multinomialResamplingAlgorithm,channelOrderCoefficientsMeans[iTrueChannelOrder+2],channelOrderCoefficientsVariances[iTrueChannelOrder+2],ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    algorithms.push_back(new CMEBasedAlgorithm("CME based algorithm (RLS no forget)",*alphabet,L,N,lastSymbolVectorInstant,noForgetRLSchannelEstimators,preamble,preamble.cols(),symbols));

    algorithms.push_back(new TimeVaryingChannelCMEbasedAlgorithm("TimeVaryingChannelCMEbasedAlgorithm",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,preamble,preamble.cols(),symbols));

    // -----------------------------------------------
//     algorithms.push_back(new LinearFilterBasedCMEapplyingAlgorithm("Linear Filter with CME based order estimation", *alphabet,L,N,lastSymbolVectorInstant,kalmanChannelEstimators,preamble, MMSEsmallLinearDetectors, ARcoefficients[0], true));

//     algorithms.push_back(new KnownSymbolsCMEapplyingAlgorithm("KnownSymbolsCMEapplyingAlgorithm",*alphabet,L,N,lastSymbolVectorInstant,noForgetRLSchannelEstimators,preamble,symbols));
    // -----------------------------------------------

//     algorithms.push_back(new ISIS("ISIS",*alphabet,L,N,lastSymbolVectorInstant,kalmanChannelEstimators,preamble,preamble.cols(),d,nParticles,multinomialResamplingAlgorithm));

    algorithms.push_back(new USIS("USIS",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,RMMSElinearDetectors,preamble,preamble.cols(),d,nParticles,multinomialResamplingAlgorithm,channelOrderEstimator,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    algorithms.push_back(new USIS2SISAlgorithm("USIS2SIS",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,RMMSElinearDetectors,preamble,preamble.cols(),d,nParticles,multinomialResamplingAlgorithm,channelOrderEstimator,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance,USISmaximumProbabilityCriterion));

    algorithms.push_back(new MLSDmAlgorithm("MLSD-m (RLS)",*alphabet,L,N,lastSymbolVectorInstant,RLSchannelEstimators,preamble,preamble.cols(),d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    algorithms.push_back(new MLSDmAlgorithm("MLSD-m (Kalman)",*alphabet,L,N,lastSymbolVectorInstant,kalmanChannelEstimators,preamble,preamble.cols(),d,nParticles,bestParticlesResamplingAlgorithm,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

    algorithms.push_back(new PSPAlgorithm("PSPAlgorithm",*alphabet,L,N,lastSymbolVectorInstant,m,kalmanEstimator,preamble,d,lastSymbolVectorInstant+d,ARcoefficients[0],nSurvivors));

    algorithms.push_back(new ViterbiAlgorithm("Viterbi",*alphabet,L,N,lastSymbolVectorInstant,*(dynamic_cast<StillMemoryMIMOChannel *> (channel)),preamble,d));
}

void TesisOrdenCanalDesconocidoSystem::BeforeEndingFrame(int iFrame)
{
    ChannelOrderEstimationSystem::BeforeEndingFrame(iFrame);
    Util::ScalarToOctaveFileStream(nSurvivors,"nSurvivors",f);
    Util::ScalarToOctaveFileStream(forgettingFactor,"forgettingFactor",f);
    Util::ScalarToOctaveFileStream(forgettingFactorDetector,"forgettingFactorDetector",f);
}

