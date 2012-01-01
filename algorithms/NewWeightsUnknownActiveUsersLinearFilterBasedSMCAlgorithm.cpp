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
#include "NewWeightsUnknownActiveUsersLinearFilterBasedSMCAlgorithm.h"

// it defines the smaller double
#include <float.h>

// #define IMPORT_REAL_DATA

// #define DEBUG

#ifdef IMPORT_REAL_DATA
    extern MIMOChannel *realChannel;
    extern MatrixXd *realSymbols;
    extern Noise *realNoise;
#endif

NewWeightsUnknownActiveUsersLinearFilterBasedSMCAlgorithm::NewWeightsUnknownActiveUsersLinearFilterBasedSMCAlgorithm(string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, LinearDetector *linearDetector, MatrixXd preamble, uint smoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm, const MatrixXd& channelMatrixMean, const MatrixXd& channelMatrixVariances, const std::vector<UsersActivityDistribution> usersActivityPdfs): SMCAlgorithm(name, alphabet, L, Nr, N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances),_linearDetector(linearDetector->clone()),_usersActivityPdfs(usersActivityPdfs)
{
    _randomParticlesInitilization = true; 
}


NewWeightsUnknownActiveUsersLinearFilterBasedSMCAlgorithm::~NewWeightsUnknownActiveUsersLinearFilterBasedSMCAlgorithm()
{
    delete _linearDetector;
}


void NewWeightsUnknownActiveUsersLinearFilterBasedSMCAlgorithm::initializeParticles()
{
      ChannelMatrixEstimator *channelMatrixEstimatorClone;

    // memory is reserved
    for(uint iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
    {
        channelMatrixEstimatorClone = _channelEstimator->clone();
        
        if(_randomParticlesInitilization)
            channelMatrixEstimatorClone->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::randnMatrix(_channelMean,_channelCovariance),rowwise,_Nr));
        
        _particleFilter->addParticle(new ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers(1.0/(double)_particleFilter->capacity(),_nInputs,_iLastSymbolVectorToBeDetected,channelMatrixEstimatorClone,_linearDetector->clone()));

        // if there is preamble...
        if(_preamble.cols()!=0)
            _particleFilter->getParticle(iParticle)->setSymbolVectors(0,_preamble.cols(),_preamble);
    }    
}

void NewWeightsUnknownActiveUsersLinearFilterBasedSMCAlgorithm::process(const MatrixXd& observations, vector< double > noiseVariances)
{
    uint iParticle,iExtendedAlphabet,iSampled;
    double s2q,sumProb,likelihood;
    MatrixXd forWeightUpdateNeededSymbols(_nInputs,_channelOrder+_d);
    
    uint iOutput,iInput;
	
	// used to store the product of the normalization constants of the proposal for the different users (needed to compute the weight)
	double normConstantsProduct;
    
	Alphabet extendedAlphabet = _alphabet.buildNewAlphabetByAddingSymbol(0.0);
	
    MatrixXd channelMatrixSample,symbolProb(_nInputs,extendedAlphabet.length());
	
	// it stores the probability of a soft estimation given each one of the symbols of the extended alphabet
	VectorXd probSoftEstGivenSymbol(extendedAlphabet.length());
	
	// it stores the product of the probability of a soft estimation given the SAMPLED symbol for all the users
	double probSoftEstGivenSampledSymbolsProduct;
	
    VectorXd sampledVector(_nInputs);    
    MatrixXd noiseCovariance = MatrixXd::Zero(_nOutputs,_nOutputs);

    for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
    {
#ifdef DEBUG
            cout << "------ iObservationToBeProcessed = " << iObservationToBeProcessed << " -----------" << endl;
#endif    
        // noise covariance needs to be constructed
        for(iOutput=0;iOutput<_nOutputs;iOutput++)
            noiseCovariance(iOutput,iOutput) = noiseVariances[iObservationToBeProcessed];

        for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
        {
#ifdef DEBUG
			cout << "iParticle = " << iParticle << endl;
#endif
            ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers *processedParticle = dynamic_cast<ParticleWithChannelEstimationAndLinearDetectionAndActiveUsers *> (_particleFilter->getParticle(iParticle));            
            
            channelMatrixSample = (dynamic_cast<KalmanEstimator *> (processedParticle->getChannelMatrixEstimator(_estimatorIndex)))->sampleFromPredictive();            
          
            // the sampled channel matrix is used to obtain soft estimations of the transmitted symbols
            VectorXd softEstimations =  processedParticle->getLinearDetector(_estimatorIndex)->detect(observations.col(iObservationToBeProcessed),channelMatrixSample,noiseCovariance);

			normConstantsProduct = 1.0;
			probSoftEstGivenSampledSymbolsProduct = 1.0;
			
            // sampling
            for(iInput=0;iInput<_nInputs;iInput++)
            {
                s2q = processedParticle->getLinearDetector(_estimatorIndex)->nthSymbolVariance(iInput,noiseVariances[iObservationToBeProcessed]);
                   
                sumProb = 0.0;
				
                // the probability for each posible symbol in the EXTENDED alphabet is computed
                for(iExtendedAlphabet=0;iExtendedAlphabet<extendedAlphabet.length();iExtendedAlphabet++)
                {
                    probSoftEstGivenSymbol(iExtendedAlphabet) = StatUtil::normalPdf(softEstimations(iInput),processedParticle->getLinearDetector(_estimatorIndex)->nthSymbolGain(iInput)*extendedAlphabet[iExtendedAlphabet],s2q);
					
					// the probability of a gaussian should never be zero (assuming the variance, "s2q", is not)
					if(probSoftEstGivenSymbol(iExtendedAlphabet)==0.0)
						probSoftEstGivenSymbol(iExtendedAlphabet) = DBL_MIN;
					
					symbolProb(iInput,iExtendedAlphabet) = probSoftEstGivenSymbol(iExtendedAlphabet);
					
					// the first time instant...
					if(iObservationToBeProcessed==_startDetectionTime)
						symbolProb(iInput,iExtendedAlphabet)*= probSymbol(extendedAlphabet[iExtendedAlphabet],_usersActivityPdfs[iInput]);
					// after first time instant
					else
						symbolProb(iInput,iExtendedAlphabet)*= probSymbolGivenPreviousActivity(extendedAlphabet[iExtendedAlphabet],processedParticle->getUserActivity(iInput,iObservationToBeProcessed-1),_usersActivityPdfs[iInput]);

                    // the computed pdf is accumulated for normalizing purposes
                    sumProb += symbolProb(iInput,iExtendedAlphabet);
                }
                
                // the product of the normalization constants for all the users is computed (sumProb is the last normalization constant computed)
                normConstantsProduct *= sumProb;

                if(sumProb!=0)
                    for(iExtendedAlphabet=0;iExtendedAlphabet<extendedAlphabet.length();iExtendedAlphabet++)
                        symbolProb(iInput,iExtendedAlphabet) /= sumProb;
                else
                {
                    for(iExtendedAlphabet=0;iExtendedAlphabet<extendedAlphabet.length();iExtendedAlphabet++)
                        symbolProb(iInput,iExtendedAlphabet) = 1.0/double(extendedAlphabet.length());
                }

                iSampled = StatUtil::discrete_rnd(symbolProb.row(iInput));
                sampledVector(iInput) = extendedAlphabet[iSampled];

				// the product of the soft estimation of the symbols given the SAMPLED symbols is updated
				probSoftEstGivenSampledSymbolsProduct *= probSoftEstGivenSymbol(iSampled);
				
				// we obtain whether the user is active or not according to the sampled symbol
				processedParticle->setUserActivity(iInput,iObservationToBeProcessed,Util::isUserActive(sampledVector(iInput)));
            } // for(iInput=0;iInput<_nInputs;iInput++)

            // sampled symbol vector is stored for the corresponding particle
            processedParticle->setSymbolVector(iObservationToBeProcessed,sampledVector);

            likelihood = StatUtil::normalPdf(observations.col(iObservationToBeProcessed),channelMatrixSample*sampledVector,noiseVariances[iObservationToBeProcessed]);

			// ********************************************* new weight computation ***********************************
			uint nCombinations = (uint) pow((double)(extendedAlphabet.length()),(double)_nInputs);

			// a likelihood is computed for every possible symbol vector
			vector<double> testedCombination(_nInputs);
			
			double combinationProb = 0.0;
			
			for(uint iTestedCombination=0;iTestedCombination<nCombinations;iTestedCombination++)
			{
				extendedAlphabet.int2symbolsArray(iTestedCombination,testedCombination);
				
// 				cout << "tested combination" << endl << testedCombination << endl;
				
				double exponent = 0.0;
				double testedCombinationPriorProb = 1.0;
				
				for(iInput=0;iInput<_nInputs;iInput++)
				{
					double mu = processedParticle->getLinearDetector(_estimatorIndex)->nthSymbolGain(iInput);
					double s2q = processedParticle->getLinearDetector(_estimatorIndex)->nthSymbolVariance(iInput,noiseVariances[iObservationToBeProcessed]);
					
// 					exponent += (pow(sampledVector(iInput),2.0) - pow(testedCombination[iInput],2.0) + 2*softEstimations(iInput)*(testedCombination[iInput]-sampledVector(iInput))) / s2q;
					exponent += mu/(2*s2q)*( mu*(pow(sampledVector(iInput),2.0) - pow(testedCombination[iInput],2.0)) + 2*softEstimations(iInput)*(testedCombination[iInput]-sampledVector(iInput)) );
					
					if(iObservationToBeProcessed==_startDetectionTime)
						testedCombinationPriorProb *= probSymbol(testedCombination[iInput],_usersActivityPdfs[iInput]);
					else
						testedCombinationPriorProb *= probSymbolGivenPreviousActivity(testedCombination[iInput],processedParticle->getUserActivity(iInput,iObservationToBeProcessed-1),_usersActivityPdfs[iInput]);
				}
				
				combinationProb += exp(exponent)*testedCombinationPriorProb;
			}

			// *****************************************************************************************************

            // the weight is updated...
// 			processedParticle->setWeight(likelihood*normConstantsProduct/probSoftEstGivenSampledSymbolsProduct *processedParticle->getWeight());
			processedParticle->setWeight(likelihood*combinationProb*processedParticle->getWeight());

			// ...the estimation of the channel matrix is updated
			processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),sampledVector,noiseVariances[iObservationToBeProcessed]);

			// ...and the channel matrix coefficients (NOT the channel matrix) stored
            processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(_estimatorIndex)->lastEstimatedChannelCoefficientsMatrix());
			
#ifdef DEBUG
			cout << "sampled symbols vector" << endl << sampledVector << endl << "updated weight = " << processedParticle->getWeight() << endl;
			cout << "the updated channel is" << endl << processedParticle->getChannelMatrix(_estimatorIndex,iObservationToBeProcessed) << endl;
#endif

        } // for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)

        _particleFilter->normalizeWeights();

		bool resamplingOccurred = false;
		
        // if it's not the last time instant
        if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
            resamplingOccurred = _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);
		
#ifdef DEBUG
		cout << "--------------- end of observation " << iObservationToBeProcessed << " ------------" << endl;
		if(!resamplingOccurred)
			cout << "NO resampling occurred: the weights after normalization:" << endl << _particleFilter->getWeightsVector() << endl;
		
		if(iObservationToBeProcessed>100)
			getchar();
#endif

    } // for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

double NewWeightsUnknownActiveUsersLinearFilterBasedSMCAlgorithm::probSymbol(tSymbol symbol,UsersActivityDistribution userActivityDistribution) const
{
	if(Util::isUserActive(symbol))
		return 1.0/double(_alphabet.length())*userActivityDistribution.probApriori(true);
	else
		return userActivityDistribution.probApriori(false);
}

double NewWeightsUnknownActiveUsersLinearFilterBasedSMCAlgorithm::probSymbolGivenPreviousActivity(tSymbol symbol,bool previousActivity,UsersActivityDistribution userActivityDistribution) const
{
	if(Util::isUserActive(symbol))
		return userActivityDistribution.probXgivenY(true,previousActivity)*1.0/double(_alphabet.length());
	else
		return userActivityDistribution.probXgivenY(false,previousActivity);
}