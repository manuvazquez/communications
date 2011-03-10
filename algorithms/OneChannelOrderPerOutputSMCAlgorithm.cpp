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

#include "OneChannelOrderPerOutputSMCAlgorithm.h"

#include <ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp.h>
#include <ChannelOrderEstimator.h>

// #define DEBUG
#define IMPORT_CHANNEL_ORDER

#ifdef IMPORT_CHANNEL_ORDER
	extern int realChannelOrder;
#endif

OneChannelOrderPerOutputSMCAlgorithm::OneChannelOrderPerOutputSMCAlgorithm(string name, Alphabet alphabet, int L, int Nr, int N, int iLastSymbolVectorToBeDetected, std::vector< ChannelMatrixEstimator* > channelEstimators, MatrixXd preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm)
:UnknownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, preamble, iFirstObservation)
,_resamplingAlgorithm(resamplingAlgorithm),_smoothingLag(smoothingLag),_randomParticlesInitilization(false)
,_channelMatrixEstimators(_nOutputs,std::vector<ChannelMatrixEstimator*>(_candidateOrders.size())),_particleFilter(new ParticleFilter(nParticles))
,_startDetectionTime(_preamble.cols())
{
	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
		if(channelEstimators[iChannelOrder]->rows()!=1)
			throw RuntimeException("OneChannelOrderPerOutputSMCAlgorithm::OneChannelOrderPerOutputSMCAlgorithm: one of the channel matrix estimators is not meant for a SISO channel.");

	for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
		for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			_channelMatrixEstimators[iOutput][iChannelOrder] = channelEstimators[iChannelOrder]->clone();

	// FIXME: this is a little bit confusing: _channelOrderAPPs[0] was correctly initialized in an upper class for the case when there is only one channel order
	_channelOrderAPPs = std::vector<MatrixXd>(_nOutputs,_channelOrderAPPs[0]);
}

OneChannelOrderPerOutputSMCAlgorithm::~OneChannelOrderPerOutputSMCAlgorithm()
{
  delete _particleFilter;
  
  for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
	  delete _channelMatrixEstimators[iOutput][iChannelOrder];
}

vector<MatrixXd> OneChannelOrderPerOutputSMCAlgorithm::getEstimatedChannelMatrices()
{
//   std::vector<MatrixXd> emptyVector;

    vector<MatrixXd> channelMatrices;
    channelMatrices.reserve(_iLastSymbolVectorToBeDetected-_preamble.cols());

    // best particle is chosen
	ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp * particle = dynamic_cast<ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp *> (_particleFilter->getBestParticle());

    for(int i=_preamble.cols();i<_iLastSymbolVectorToBeDetected;i++)
	{
	  MatrixXd channelMatrix = MatrixXd::Zero(_nOutputs,_nInputsXmaxChannelOrder);
	  for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
	  {
// 		cout << "best channel order = " << particle->iMaxChannelOrderAPP(iOutput) << endl;
		MatrixXd outputChannelMatrix = particle->getChannelMatrix(iOutput,particle->iMaxChannelOrderAPP(iOutput),i);
		
		if(outputChannelMatrix.rows()>1)
		  throw RuntimeException("OneChannelOrderPerOutputSMCAlgorithm::getEstimatedChannelMatrices: estimated channel for one output is bigger than 1.");
		
		channelMatrix.block(iOutput,_nInputsXmaxChannelOrder-outputChannelMatrix.cols(),1,outputChannelMatrix.cols()) = outputChannelMatrix;
		
// 		cout << "channelMatrix = " << endl << channelMatrix << endl;
// 		cout << "est" << endl << outputChannelMatrix << endl;
	  }
	  channelMatrices.push_back(channelMatrix);
	}

    return channelMatrices;  
//   return emptyVector;
}

MatrixXd OneChannelOrderPerOutputSMCAlgorithm::getDetectedSymbolVectors()
{
  return particleFilter()->getBestParticle()->getSymbolVectors(_preamble.cols(),_iLastSymbolVectorToBeDetected-1);
}

void OneChannelOrderPerOutputSMCAlgorithm::run(MatrixXd observations, std::vector<double> noiseVariances, MatrixXd trainingSequence)
{
	std::vector<std::vector<MatrixXd> >  estimatedChannelMatrices = processTrainingSequence(observations,noiseVariances,trainingSequence);

	_startDetectionTime += trainingSequence.cols();

	initializeParticles();

	ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp * initialParticle = dynamic_cast<ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp *> (_particleFilter->getParticle(0));

	for(int i=_preamble.cols();i<_preamble.cols()+trainingSequence.cols();i++)
		for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
			for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			initialParticle->setChannelMatrix(iOutput,iChannelOrder,i,estimatedChannelMatrices[iOutput][iChannelOrder]);

	initialParticle->setSymbolVectors(_preamble.cols(),_preamble.cols()+trainingSequence.cols(),trainingSequence);


#ifdef DEBUG
	for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
		for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			cout << "channel order probability " << _candidateOrders[iChannelOrder] << " for output " << iOutput << " at " <<  _preamble.cols()+trainingSequence.cols()-1 << ": " << _channelOrderAPPs[iOutput](iChannelOrder,_preamble.cols()+trainingSequence.cols()-1) << endl;

	cout << "before process" << endl;
#endif

	process(observations,noiseVariances);

#ifdef DEBUG
	cout << "finished run!!" << endl;
	//   getchar();
#endif
}

void OneChannelOrderPerOutputSMCAlgorithm::run(MatrixXd observations, std::vector<double> noiseVariances)
{
	throw RuntimeException("OneChannelOrderPerOutputSMCAlgorithm::run: this algorithm is not implemented to run without training sequence.");
}

std::vector<std::vector<MatrixXd> > OneChannelOrderPerOutputSMCAlgorithm::processTrainingSequence(const MatrixXd &observations, const std::vector<double> &noiseVariances, const MatrixXd &trainingSequence)
{
    MatrixXd sequenceToProcess(trainingSequence.rows(),_preamble.cols()+trainingSequence.cols());
    sequenceToProcess << _preamble,trainingSequence;

	// if there is no enough observations to process the training sequence...
    if(observations.cols() < (_preamble.cols()+trainingSequence.cols()))
        throw RuntimeException("OneChannelOrderPerOutputSMCAlgorithm::run: not enough number of observations to process the training sequence.");

#ifdef DEBUG
	cout << _channelOrderAPPs.size() << endl;
#endif
	
  for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
	for(uint iChannelOrder=0;iChannelOrder<static_cast<uint>(_candidateOrders.size());iChannelOrder++)
	  // at the beginning, all the channel orders have the same probability for all the channels
	  _channelOrderAPPs[iOutput](iChannelOrder,_preamble.cols()-1) = 1.0/double(_candidateOrders.size());

    double normConst;
    vector<double> unnormalizedChannelOrderAPPs(_candidateOrders.size());

#ifdef DEBUG
  cout << "_nOutputs = " << _nOutputs << endl;
#endif

  std::vector<std::vector<MatrixXd> > estimatedChannelMatrices = std::vector<std::vector<MatrixXd> >(_nOutputs,std::vector<MatrixXd>(_candidateOrders.size()));

  for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
  {
#ifdef DEBUG
  cout << "===================================iOutput = " << iOutput << " =======================" << endl;
#endif
	for(int i=_preamble.cols();i<sequenceToProcess.cols();i++)
	{
#ifdef DEBUG
	cout << " ---------------------------- time instant " << i << endl;
#endif
	  normConst = 0.0;
	  for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
	  {
		
#ifdef DEBUG
		  cout << "iChannelOrder = " << iChannelOrder << " " << _channelMatrixEstimators[iOutput][iChannelOrder]->likelihood(observations.block(iOutput,i,1,1),sequenceToProcess.block(0,i-_candidateOrders[iChannelOrder]+1,_nInputs,_candidateOrders[iChannelOrder]),noiseVariances[i]) << endl;
#endif
		  // unnormalized channel order APP
		  unnormalizedChannelOrderAPPs[iChannelOrder] = _channelOrderAPPs[iOutput](iChannelOrder,i-1)*_channelMatrixEstimators[iOutput][iChannelOrder]->likelihood(observations.block(iOutput,i,1,1),sequenceToProcess.block(0,i-_candidateOrders[iChannelOrder]+1,_nInputs,_candidateOrders[iChannelOrder]),noiseVariances[i]);
		  normConst += unnormalizedChannelOrderAPPs[iChannelOrder];

		  _channelMatrixEstimators[iOutput][iChannelOrder]->nextMatrix(observations.block(iOutput,i,1,1),sequenceToProcess.block(0,i-_candidateOrders[iChannelOrder]+1,_nInputs,_candidateOrders[iChannelOrder]),noiseVariances[i]);
		  estimatedChannelMatrices[iOutput][iChannelOrder] = _channelMatrixEstimators[iOutput][iChannelOrder]->lastEstimatedChannelMatrix();
	  }

	  if(normConst!=0.0)
		  for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			  _channelOrderAPPs[iOutput](iChannelOrder,i) = unnormalizedChannelOrderAPPs[iChannelOrder] / normConst;
	  else
		  for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			  _channelOrderAPPs[iOutput](iChannelOrder,i) = _channelOrderAPPs[iOutput](iChannelOrder,i-1);
#ifdef DEBUG
	  cout << _channelOrderAPPs[iOutput](0,i) << endl;
#endif
	}
  }
  
  return estimatedChannelMatrices;
}

void OneChannelOrderPerOutputSMCAlgorithm::initializeParticles()
{
    std::vector<std::vector<ChannelMatrixEstimator*> > channelEstimatorsClone(_nOutputs,std::vector<ChannelMatrixEstimator*>(_candidateOrders.size()));
	
  for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
	for(uint iChannelOrder=0;iChannelOrder<static_cast<uint>(_candidateOrders.size());iChannelOrder++)
	  channelEstimatorsClone[iOutput][iChannelOrder] = _channelMatrixEstimators[iOutput][iChannelOrder]->clone();

    // we begin with only one particle
	ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp *particle = new ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp(1.0,_nInputs,_iLastSymbolVectorToBeDetected+_smoothingLag,channelEstimatorsClone);

    particle->setSymbolVectors(0,_preamble.cols(),_preamble);

    // the available APP's just before the _startDetectionTime instant are copied into the particle
	for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
	  for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
	  {
#ifdef DEBUG
		cout << "iOutput = " << iOutput << " iChannelOrder = " << iChannelOrder << endl;
#endif
		particle->setChannelOrderAPP(iOutput,_channelOrderAPPs[iOutput](iChannelOrder,_startDetectionTime-1),iChannelOrder);
	  }

    _particleFilter->addParticle(particle);
}

void OneChannelOrderPerOutputSMCAlgorithm::process(const MatrixXd &observations,const vector<double> &noiseVariances)
{
    uint nSymbolVectors = (int) pow((double)_alphabet.length(),(double)_nInputs);
    vector<tSymbol> testedVector(_nInputs);
    VectorXd computedObservations(_nOutputs);
    int iCandidate,m,k,iParticle;
    uint iChannelOrder,iTestedVector;
    ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp *processedParticle;
	
	// over all the survivor candidates
    double normConst;
	
	double singleOutputLikelihood;

//     typedef struct{
//         int fromParticle;
//         MatrixXd symbolVectorsMatrix;
// 		MatrixXd channelOrderAPPs;
//         double weight;
//     }tParticleCandidate;

    tParticleCandidate *particleCandidates = new tParticleCandidate[_particleFilter->capacity()*nSymbolVectors];

    // it will contain all the symbols involved in the current observation
    MatrixXd symbolVectorsMatrix(_nInputs,_maxOrder);
	
    VectorXd symbolsVector;

	// when all the symbols involved in an observation are stacked row-wise, this indicates where the last symbol vector starts
    int iLastSymbolVectorStartWithinVector = _nInputsXmaxChannelOrder - _nInputs;

    for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected+_smoothingLag;iObservationToBeProcessed++)
    {
#ifdef DEBUG
		cout << "============= iObservationToBeProcessed = " << iObservationToBeProcessed << " ==============" << endl;
#endif
        // it keeps track of the place where a new tParticleCandidate will be stored within the array
        iCandidate = 0;

        normConst = 0.0;

        // the candidates from all the particles are generated
        for(iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)
        {
            processedParticle = dynamic_cast<ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp *> (_particleFilter->getParticle(iParticle));

            symbolVectorsMatrix.block(0,0,_nInputs,_maxOrder-1) = processedParticle->getSymbolVectors(iObservationToBeProcessed-_maxOrder+1,iObservationToBeProcessed-1);
            symbolsVector = Util::toVector(symbolVectorsMatrix,columnwise);

            for(iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)
            {
                // the corresponding testing vector is generated from the index
                _alphabet.int2symbolsArray(iTestedVector,testedVector);

                // current tested vector is copied in the m-th position
                for(k=0;k<_nInputs;k++)
                    symbolVectorsMatrix(k,_maxOrder-1) = symbolsVector(iLastSymbolVectorStartWithinVector+k) = testedVector[k];

				double vectorLikelihood = 1.0;
				
				std::vector<double> unnormalizedChannelOrderAPPs(_candidateOrders.size());
				
				particleCandidates[iCandidate].channelOrderAPPs = MatrixXd(_nOutputs,_candidateOrders.size());

				for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
				{
#ifdef DEBUG
				  cout << "......... iOutput = " << iOutput << " .............." << endl;
#endif
				  singleOutputLikelihood = 0.0;
				
				  for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
				  {
					  m = _candidateOrders[iChannelOrder];

					  MatrixXd involvedSymbolVectors = symbolVectorsMatrix.block(0,_maxOrder-m,_nInputs,m);

					  unnormalizedChannelOrderAPPs[iChannelOrder] = processedParticle->getChannelOrderAPP(iOutput,iChannelOrder)*
									processedParticle->getChannelMatrixEstimator(iOutput,iChannelOrder)->likelihood(observations.block(iOutput,iObservationToBeProcessed,1,1),involvedSymbolVectors,noiseVariances[iObservationToBeProcessed]);


					  singleOutputLikelihood += unnormalizedChannelOrderAPPs[iChannelOrder];
#ifdef DEBUG
					  cout << " //////// channelOrder = " << _candidateOrders[iChannelOrder] << " \\\\\\\\\\\\" << endl;
					  cout << processedParticle->getChannelMatrixEstimator(iOutput,iChannelOrder)->lastEstimatedChannelMatrix() << endl;
					  getchar();
#endif
				  }
				  
				  // normalization
				  for(iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
					particleCandidates[iCandidate].channelOrderAPPs(iOutput,iChannelOrder) = unnormalizedChannelOrderAPPs[iChannelOrder]/singleOutputLikelihood;
				  
				  vectorLikelihood *= singleOutputLikelihood;
				} // for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)

                // if the likelihood is zero, we don't generate a candidate for this particle and this symbol vector
                if(vectorLikelihood==0.0)
                    continue;

                particleCandidates[iCandidate].fromParticle = iParticle;
                particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;
                particleCandidates[iCandidate].weight = processedParticle->getWeight()*vectorLikelihood;
				
				// normalization constant is computed taking advantage of the loop
                normConst += particleCandidates[iCandidate].weight;

                iCandidate++;
            } // for(uint iTestedVector=0;iTestedVector<nSymbolVectors;iTestedVector++)

        } // for(int iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)

        // if none of the candidates was valid
        if(iCandidate==0)
        {
            VectorXd uniformDistribution(_alphabet.length());
            for(int iAlphabet=0;iAlphabet<_alphabet.length();iAlphabet++)
                uniformDistribution(iAlphabet) = 1.0/_alphabet.length();

            for(iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)
            {
                processedParticle = dynamic_cast<ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp *> (_particleFilter->getParticle(iParticle));
                
                symbolVectorsMatrix.block(0,0,_nInputs,_maxOrder-1) = processedParticle->getSymbolVectors(iObservationToBeProcessed-_maxOrder+1,iObservationToBeProcessed-1);
                
				// the transmitted symbols are sampled from a uniform distribution
                for(k=0;k<_nInputs;k++)
                    symbolVectorsMatrix(k,_maxOrder-1) = _alphabet[StatUtil::discrete_rnd(uniformDistribution)];
                
				particleCandidates[iParticle].fromParticle = iParticle;
                particleCandidates[iCandidate].symbolVectorsMatrix = symbolVectorsMatrix;

                particleCandidates[iCandidate].weight = processedParticle->getWeight();
				
				for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
				  for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
					particleCandidates[iCandidate].channelOrderAPPs(iOutput,iChannelOrder) = processedParticle->getChannelOrderAPP(iOutput,iChannelOrder);
				
				// normalization constant is computed taking advantage of the loop
                normConst += particleCandidates[iCandidate].weight;
            }
            iCandidate = _particleFilter->nParticles();
        }

        // a vector of size the number of generated candidates is declared...
        VectorXd weights(iCandidate);

        // ...to store their weights
        for(int i=0;i<iCandidate;i++)
            weights(i) = particleCandidates[i].weight/normConst;

		// an overall number of survivors is considered
		// the candidates that are going to give rise to particles are selected
        vector<int> indexesSelectedCandidates = _resamplingAlgorithm->obtainIndexes(_particleFilter->capacity(),weights);

		// ------------------------------ fixed number of survivors per state -----------------------------------------------------------
		
// 		vector<int> indexesSelectedCandidates;
// 		
// 		uint nSurvivors = _particleFilter->capacity()/int(pow(double(_alphabet.length()),double(_nInputs*(realChannelOrder-1))));
// 		if(_particleFilter->capacity() % int(pow(double(_alphabet.length()),double(_nInputs*(realChannelOrder-1)))) !=0)
// 			throw RuntimeException("OneChannelOrderPerOutputSMCAlgorithm:process: the number of computed survivors is not integer.");
// 		
// 		std::vector<std::vector<bool> > stateMasks = imposeFixedNumberOfSurvivorsPerState(particleCandidates,iCandidate);
// 		for(uint i=0;i<stateMasks.size();i++)
// 		{
// 			std::vector<int> thisStateIndexes = _resamplingAlgorithm->obtainIndexes(nSurvivors,weights,stateMasks[i]);
// 			indexesSelectedCandidates.insert(indexesSelectedCandidates.end(),thisStateIndexes.begin(),thisStateIndexes.end());
// 		}

		// ------------------------------ fixed number of survivors per state -----------------------------------------------------------

        // every survivor candidate is associated with an old particle
        vector<int> indexesParticles(indexesSelectedCandidates.size());
        for(uint i=0;i<indexesSelectedCandidates.size();i++)
            indexesParticles[i] = particleCandidates[indexesSelectedCandidates[i]].fromParticle;

        // the chosen particles are kept without modification (yet)
        _particleFilter->keepParticles(indexesParticles);

        // every surviving particle is modified according to what it says its corresponding candidate
        for(int iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)
        {
            processedParticle = dynamic_cast<ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp *> (_particleFilter->getParticle(iParticle));

            // sampled symbols are copied into the corresponding particle
            processedParticle->setSymbolVector(iObservationToBeProcessed,particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.col(_maxOrder-1));

			for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
			{
			  for(uint iChannelOrder=0;iChannelOrder<processedParticle->nChannelMatrixEstimators(iOutput);iChannelOrder++)
			  {
				  // channel matrix is estimated by means of the particle channel estimator
				  processedParticle->setChannelMatrix(iOutput,iChannelOrder,iObservationToBeProcessed,processedParticle->
					getChannelMatrixEstimator(iOutput,iChannelOrder)->
					  nextMatrix(observations.block(iOutput,iObservationToBeProcessed,1,1),particleCandidates[indexesSelectedCandidates[iParticle]].symbolVectorsMatrix.block(0,_maxOrder-_candidateOrders[iChannelOrder],_nInputs,_candidateOrders[iChannelOrder]),noiseVariances[iObservationToBeProcessed]));

				  // notice that the likelihood of a candidate is the sum of the likelihoods for the different channel orders, i.e, the normalization constant
				  processedParticle->setChannelOrderAPP(iOutput,particleCandidates[indexesSelectedCandidates[iParticle]].channelOrderAPPs(iOutput,iChannelOrder),iChannelOrder);
			  }
			}

            processedParticle->setWeight(particleCandidates[indexesSelectedCandidates[iParticle]].weight);

        } // for(int iParticle=0;iParticle<_particleFilter->nParticles();iParticle++)

        _particleFilter->normalizeWeights();

        processedParticle = dynamic_cast<ParticleWithMultipleChannelsEstimationAndMultipleChannelOrderApp *> (_particleFilter->getBestParticle());

		for (uint iOutput=0;iOutput<static_cast<uint>(_nOutputs);iOutput++)
		  for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
			_channelOrderAPPs[iOutput](iChannelOrder,iObservationToBeProcessed) = processedParticle->getChannelOrderAPP(iOutput,iChannelOrder);

    } // for(int iObservationToBeProcessed=_startDetectionTime;iObservationToBeProcessed<_iLastSymbolVectorToBeDetected+_d;iObservationToBeProcessed++)

    delete[] particleCandidates;
}

std::vector<std::vector<bool> > OneChannelOrderPerOutputSMCAlgorithm::imposeFixedNumberOfSurvivorsPerState(const tParticleCandidate *particleCandidates,uint nCandidates)
{
	int nStates = pow(double(_alphabet.length()),double(_nInputs*(realChannelOrder-1)));
	
	std::vector<std::vector<bool> > statesMasks(nStates,std::vector<bool>(nCandidates,false));
	
	for(int iState=0;iState<nStates;iState++)
	{
		MatrixXd state = _alphabet.int2eigenMatrix(iState,_nInputs,realChannelOrder-1);
		
		for(uint iCandidate=0;iCandidate<nCandidates;iCandidate++)
		{
			if(particleCandidates[iCandidate].symbolVectorsMatrix.block(0,particleCandidates[iCandidate].symbolVectorsMatrix.cols()-realChannelOrder+1,_nInputs,realChannelOrder-1) == state)
			{
				statesMasks[iState][iCandidate] = true;
			}
		}
	}
	
	return statesMasks;
}