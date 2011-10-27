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

#ifndef ONECHANNELORDERPEROUTPUTSMCALGORITHM_H
#define ONECHANNELORDERPEROUTPUTSMCALGORITHM_H

#include <UnknownChannelOrderAlgorithm.h>

#include <ResamplingAlgorithm.h>

class OneChannelOrderPerOutputSMCAlgorithm : public UnknownChannelOrderAlgorithm
{
protected:
	
    typedef struct{
        int fromParticle;
        MatrixXd symbolVectorsMatrix;
		MatrixXd channelOrderAPPs;
        double weight;
    }tParticleCandidate;

    ResamplingAlgorithm *_resamplingAlgorithm;
    uint _smoothingLag;
    
    // indicates whether the particles will be initialized randomly or all the same
    bool _randomParticlesInitilization;

	// _channelMatrixEstimators[i][j] is the the channel matrix estimator for the j-th channel order on the i-th output
	std::vector<std::vector<ChannelMatrixEstimator*> > _channelMatrixEstimators;
	
	ParticleFilter *_particleFilter;
	uint _startDetectionTime;

	/*!
	  It returns a pointer to the particle filter used by the algorithm
	  \return a pointer to a particle filter
	*/
	virtual ParticleFilter* particleFilter() { return _particleFilter;}

	virtual void initializeParticles();

	virtual void process(const MatrixXd &observations,const vector<double> &noiseVariances);

	std::vector<std::vector<MatrixXd> > processTrainingSequence(const MatrixXd &observations, const std::vector<double> &noiseVariances, const MatrixXd &trainingSequence);
	
	std::vector<std::vector<bool> > imposeFixedNumberOfSurvivorsPerState(const tParticleCandidate *particleCandidates, uint nCandidates);

public:
  
	OneChannelOrderPerOutputSMCAlgorithm(string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected,vector<ChannelMatrixEstimator *> channelEstimators,MatrixXd preamble,uint iFirstObservation,uint smoothingLag,uint nParticles,ResamplingAlgorithm *resamplingAlgorithm);
	~OneChannelOrderPerOutputSMCAlgorithm();

	virtual bool estimatesOneChannelOrderPerOutput() const { return true;}
	virtual bool estimatesOneSingleChannelOrder() const { return false;}

    virtual vector<MatrixXd> getEstimatedChannelMatrices();
    virtual MatrixXd getDetectedSymbolVectors();
    virtual void run(MatrixXd observations, std::vector< double> noiseVariances, MatrixXd trainingSequence);
    virtual void run(MatrixXd observations, std::vector< double> noiseVariances);
};

#endif // ONECHANNELORDERPEROUTPUTSMCALGORITHM_H
