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

OneChannelOrderPerOutputSMCAlgorithm::OneChannelOrderPerOutputSMCAlgorithm(string name, Alphabet alphabet, int L, int Nr, int N, int iLastSymbolVectorToBeDetected, std::vector< ChannelMatrixEstimator* > channelEstimators, MatrixXd preamble, int iFirstObservation, int smoothingLag, int nParticles, ResamplingAlgorithm* resamplingAlgorithm)
:UnknownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, preamble, iFirstObservation)
,_resamplingAlgorithm(resamplingAlgorithm),_smoothingLag(smoothingLag),_randomParticlesInitilization(false)
// ,_channelMatrixEstimators(_L,std::<ChannelMatrixEstimator*>
{
//   for (uint iOutput=0;iOutput<_L;iOutput++)
	
}

std::vector< MatrixXd, std::allocator< MatrixXd > > OneChannelOrderPerOutputSMCAlgorithm::getEstimatedChannelMatrices()
{

}

MatrixXd OneChannelOrderPerOutputSMCAlgorithm::getDetectedSymbolVectors()
{
  return particleFilter()->getBestParticle()->getSymbolVectors(_preamble.cols(),_iLastSymbolVectorToBeDetected-1);
}

void OneChannelOrderPerOutputSMCAlgorithm::run(MatrixXd observations, std::vector< double, std::allocator< double > > noiseVariances, MatrixXd trainingSequence)
{

}

void OneChannelOrderPerOutputSMCAlgorithm::run(MatrixXd observations, std::vector< double, std::allocator< double > > noiseVariances)
{
	throw RuntimeException("OneChannelOrderPerOutputSMCAlgorithm::run: this algorithm is not implemented to run without training sequence.");
}

