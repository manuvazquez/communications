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
    ResamplingAlgorithm *_resamplingAlgorithm;
    int _smoothingLag;
//     int _startDetectionTime;

// 	// mean and variance which will server to initialize
//     double _channelUniqueMean, _channelUniqueVariance;
// 
// 	// 
//     vector<MatrixXd> _channelMeanVectors;
//     vector<MatrixXd> _channelCovariances;
    
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

//     virtual void process(const MatrixXd &observations,vector<double> noiseVariances) = 0;
	virtual void process(const MatrixXd &observations,const vector<double> &noiseVariances);

	void processTrainingSequence(const MatrixXd &observations, const std::vector<double> &noiseVariances, const MatrixXd &trainingSequence);

//     virtual int iBestChannelOrder(int iBestParticle) = 0;

//     virtual void beforeInitializingParticles(const MatrixXd &observations,vector<double> &noiseVariances,const MatrixXd &trainingSequence) {}
//     virtual void updateParticleChannelOrderEstimators(Particle *particle,const MatrixXd &observations,const std::vector<std::vector<MatrixXd> > &channelMatrices,vector<double> &noiseVariances,const MatrixXd &sequenceToProcess) {}
public:
  
	OneChannelOrderPerOutputSMCAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected,vector<ChannelMatrixEstimator *> channelEstimators,MatrixXd preamble,int iFirstObservation,int smoothingLag,int nParticles,ResamplingAlgorithm *resamplingAlgorithm);
	~OneChannelOrderPerOutputSMCAlgorithm();

	virtual bool estimatesOneChannelOrderPerOutput() const { return true;}
	virtual bool estimatesOneSingleChannelOrder() const { return false;}

    virtual vector<MatrixXd> getEstimatedChannelMatrices();
    virtual MatrixXd getDetectedSymbolVectors();
    virtual void run(MatrixXd observations, std::vector< double> noiseVariances, MatrixXd trainingSequence);
    virtual void run(MatrixXd observations, std::vector< double> noiseVariances);
};

#endif // ONECHANNELORDERPEROUTPUTSMCALGORITHM_H
