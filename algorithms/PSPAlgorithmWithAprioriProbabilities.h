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

#ifndef PSPALGORITHMWITHAPRIORIPROBABILITIES_H
#define PSPALGORITHMWITHAPRIORIPROBABILITIES_H

#include <PSPAlgorithm.h>

#include <UsersActivityDistribution.h>
#include <KalmanEstimator.h>

class PSPAlgorithmWithAprioriProbabilities : public PSPAlgorithm
{

protected:
  const std::vector<UsersActivityDistribution> _usersActivityPdfs;
  Alphabet _extendedAlphabet;

  virtual void deployState(uint iState, const VectorXd& observations, double noiseVariance);

public:

  PSPAlgorithmWithAprioriProbabilities(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, uint firstSymbolVectorDetectedAt, uint nSurvivors, const std::vector<UsersActivityDistribution> usersActivityPdfs);
  
  virtual void run(MatrixXd observations, vector< double > noiseVariances);
  virtual void run(MatrixXd observations, vector< double > noiseVariances, MatrixXd trainingSequence);
};

#endif // PSPALGORITHMWITHAPRIORIPROBABILITIES_H
