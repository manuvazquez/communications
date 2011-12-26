/*
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License version 2 as published by the Free Software Foundation.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public License
   along with this library; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/

#include "CDMAunknownActiveUsersSISopt.h"

#include <ParticleWithChannelEstimation.h>

// #define DEBUG3
// #define DEBUG4

CDMAunknownActiveUsersSISopt::CDMAunknownActiveUsersSISopt(string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm, const MatrixXd& channelMatrixMean, const MatrixXd& channelMatrixVariances,const std::vector<UsersActivityDistribution> usersActivityPdfs): SMCAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, m, channelEstimator, preamble, smoothingLag, nParticles, resamplingAlgorithm, channelMatrixMean, channelMatrixVariances),_usersActivityPdfs(usersActivityPdfs)
{    
  _randomParticlesInitilization = true;    
}

void CDMAunknownActiveUsersSISopt::process(const MatrixXd& observations, std::vector< double, std::allocator< double > > noiseVariances)
{
  // a new alphabet extended with 0 (that meaning, no symbol is transmitted)
  Alphabet extendedAlphabet = _alphabet.buildNewAlphabetByAddingSymbol(0.0);
	  
  uint nCombinations = (uint) pow((double)(extendedAlphabet.length()),(double)_nInputs);

  uint k,iParticle,iSampledVector;
  vector<tSymbol> sampledVector(_nInputs);

  // tVector containing the symbols
  VectorXd symbolsVector(_nInputs);

  // a likelihood is computed for every possible symbol vector
  vector<double> testedCombination(_nInputs),likelihoods(nCombinations);

  double likelihoodsSum;

  // for each time instant
  for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
  {
	  for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
	  {
		  ParticleWithChannelEstimation *processedParticle = dynamic_cast<ParticleWithChannelEstimation *> (_particleFilter->getParticle(iParticle));

#ifdef DEBUG
		  if(iObservationToBeProcessed!=_startDetectionTime)
		  {
			cout << "previous symbol vector" << endl;
			cout << processedParticle->getSymbolVector(iObservationToBeProcessed-1) << endl;
		  }
#endif
		  
		  for(uint iTestedCombination=0;iTestedCombination<nCombinations;iTestedCombination++)
		  {
			  // the corresponding testing vector is generated from the index
			  extendedAlphabet.int2symbolsArray(iTestedCombination,testedCombination);

			  // symbols are copied into a vector
			  for(k=0;k<_nInputs;k++)
				  symbolsVector(k) = testedCombination[k];

#ifdef DEBUG2
			  cout << "testing iTestedCombination = " << iTestedCombination << endl;
			  cout << symbolsVector << endl;
#endif
			  
			  likelihoods[iTestedCombination] = processedParticle->getChannelMatrixEstimator()->likelihood(observations.col(iObservationToBeProcessed),symbolsVector,noiseVariances[iObservationToBeProcessed]);
			  
#ifdef DEBUG4
			  cout << "variances" << endl << processedParticle->getChannelMatrixEstimator()->getVariances() << endl;
#endif

			  // the probability of these users being active and this particular symbol vector being transmitted is computed...
			  // ...either taking into account the users previous state in case this exists
			  if(iObservationToBeProcessed!=_startDetectionTime)
				  likelihoods[iTestedCombination] *= StatUtil::probSymbolsVectorGivenPreviousTimeInstantUsersActivity(
													  symbolsVector,Util::getUsersActivityFromSymbolsVector(processedParticle->getSymbolVector(iObservationToBeProcessed-1)),
													  _usersActivityPdfs,_alphabet.length());
			  // ...or not doing so when it's the first time instant
			  else
				  likelihoods[iTestedCombination] *= StatUtil::probSymbolsVector(symbolsVector,_usersActivityPdfs,_alphabet.length());

#ifdef DEBUG
			  cout << "iTestedCombination " <<  endl << symbolsVector << endl;
			  if(iObservationToBeProcessed!=_startDetectionTime)
				cout << "prob symbols = " << probSymbolsVectorGivenPreviousTimeInstantUsersActivity(symbolsVector,getUsersActivityFromSymbolsVector(processedParticle->getSymbolVector(iObservationToBeProcessed-1))) << endl;
			  else
				cout << "prob symbols (a priori) = " << probSymbolsVectorGivenPreviousTimeInstantUsersActivity(symbolsVector) << endl;
#endif
				  
		  } // for(uint iTestedCombination=0;iTestedCombination<nCombinations;iTestedCombination++)

		  likelihoodsSum = Util::sum(likelihoods);
		  vector<double> probabilities = likelihoods;
		  if(likelihoodsSum>0)
			  // probabilities are computed by normalizing the likelihoods
			  Util::normalize(probabilities);
		  else
			  // if all the likelihoods are null
			  probabilities = vector<double>(nCombinations,1.0/(double)nCombinations);

		  // one sample from the discrete distribution is taken
		  iSampledVector = StatUtil::discrete_rnd(probabilities);
		  
  //             iSampledVector = Util::max(probabilities); // <---------------------------------------------------------------------------------

		  // the above index is turned into a vector
		  extendedAlphabet.int2symbolsArray(iSampledVector,sampledVector);

#ifdef DEBUG
		  Util::print(sampledVector);
		  cout << "iSampledVector = " << iSampledVector << endl;
		  cout << "-----------------------------------------------" << endl;
		  getchar();
#endif
		  
		  // sampled symbols are copied into the corresponding particle
		  processedParticle->setSymbolVector(iObservationToBeProcessed,sampledVector);

		  // channel matrix is estimated by means of the particle channel estimator
		  processedParticle->getChannelMatrixEstimator(_estimatorIndex)->nextMatrix(observations.col(iObservationToBeProcessed),processedParticle->getSymbolVector(iObservationToBeProcessed),noiseVariances[iObservationToBeProcessed]);

		  // the last estimated channel coefficients (NOT channel matrix) are stored in the particle
		  processedParticle->setChannelMatrix(_estimatorIndex,iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(_estimatorIndex)->lastEstimatedChannelCoefficientsMatrix());
		  
#ifdef SAVE_CHANNEL_ESTIMATES_VARIANCES
		  processedParticle->setChannelEstimatesVariances(iObservationToBeProcessed,processedParticle->getChannelMatrixEstimator(_estimatorIndex)->getVariances());
#endif

		  processedParticle->setWeight(processedParticle->getWeight()* likelihoodsSum);
		  
#ifdef DEBUG3
		  cout << *processedParticle << endl;
#endif

	  } // for(iParticle=0;iParticle<_particleFilter->capacity();iParticle++)

	  _particleFilter->normalizeWeights();

	  // if it's not the last time instant
	  if(iObservationToBeProcessed<(_iLastSymbolVectorToBeDetected-1))
		  _resamplingAlgorithm->resampleWhenNecessary(_particleFilter);
	  
#ifdef DEBUG3
		  cout << "end of processing of observation " << iObservationToBeProcessed << endl;
		  cout << "pesos" << endl << _particleFilter->getWeightsVector() << endl;
		  getchar();
#endif

  } // for(uint iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected;iObservationToBeProcessed++)
}

void CDMAunknownActiveUsersSISopt::initializeParticles()
{
    ChannelMatrixEstimator *channelMatrixEstimatorClone;

    // memory is reserved
    for(uint iParticle=0;iParticle<_particleFilter->capacity();iParticle++)
    {
        channelMatrixEstimatorClone = _channelEstimator->clone();
        
        if(_randomParticlesInitilization)
            channelMatrixEstimatorClone->setFirstEstimatedChannelMatrix(Util::toMatrix(StatUtil::randnMatrix(_channelMean,_channelCovariance),rowwise,_Nr));
        
        _particleFilter->addParticle(new ParticleWithChannelEstimation(1.0/(double)_particleFilter->capacity(),_nInputs,_iLastSymbolVectorToBeDetected,channelMatrixEstimatorClone));

        // if there is preamble...
        if(_preamble.cols()!=0)
            _particleFilter->getParticle(iParticle)->setSymbolVectors(0,_preamble.cols(),_preamble);
    }
}
