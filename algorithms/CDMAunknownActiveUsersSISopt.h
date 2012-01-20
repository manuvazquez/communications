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

#ifndef CDMAUNKNOWNACTIVEUSERSSISOPT_H
#define CDMAUNKNOWNACTIVEUSERSSISOPT_H

#include <SMCAlgorithm.h>
#include <UsersActivityDistribution.h>


class CDMAunknownActiveUsersSISopt : public SMCAlgorithm
{
public:
	CDMAunknownActiveUsersSISopt(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, uint m, ChannelMatrixEstimator* channelEstimator, MatrixXd preamble, uint smoothingLag, uint nParticles, ResamplingAlgorithm* resamplingAlgorithm, const MatrixXd& channelMatrixMean, const MatrixXd& channelMatrixVariances,const std::vector<UsersActivityDistribution> usersActivityPdfs);
protected:
	const std::vector<UsersActivityDistribution> _usersActivityPdfs; /// objects describing the pdf of the users activity
  
	virtual void process(const MatrixXd& observations, std::vector< double, std::allocator< double > > noiseVariances);
    virtual void initializeParticles();
};

#endif // CDMAUNKNOWNACTIVEUSERSSISOPTWITHUSERSACTIVITYSAMPLING_H
