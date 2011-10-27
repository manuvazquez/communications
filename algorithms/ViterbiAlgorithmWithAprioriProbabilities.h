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

#ifndef VITERBIALGORITHMWITHAPRIORIPROBABILITIES_H
#define VITERBIALGORITHMWITHAPRIORIPROBABILITIES_H

#include <ViterbiAlgorithm.h>
#include <UsersActivityDistribution.h>


class ViterbiAlgorithmWithAprioriProbabilities : public ViterbiAlgorithm
{
public:
    ViterbiAlgorithmWithAprioriProbabilities(string name, Alphabet alphabet,uint L,uint Nr,uint N, uint iLastSymbolVectorToBeDetected, const StillMemoryMIMOChannel& channel,const MatrixXd &preamble,uint smoothingLag, const std::vector<UsersActivityDistribution> usersActivityPdfs);
	
	virtual void run(MatrixXd observations,vector<double> noiseVariances,uint firstSymbolVectorDetectedAt);
protected:
    const std::vector<UsersActivityDistribution> _usersActivityPdfs;
	Alphabet _extendedAlphabet;
  
  virtual void deployState(int iState, const VectorXd& observations, const MatrixXd& channelMatrix, const double noiseVariance);
};

#endif // VITERBIALGORITHMWITHAPRIORIPROBABILITIES_H
