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
#ifndef WITHCHANNELESTIMATIONPARTICLEADDON_H
#define WITHCHANNELESTIMATIONPARTICLEADDON_H

/**
It implements the channel related part for a particle

	@author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <vector>
#include <ChannelMatrixEstimator.h>

class WithChannelEstimationParticleAddon{
protected:
    std::vector<ChannelMatrixEstimator *> _channelMatrixEstimators;
    std::vector<std::vector<MatrixXd> > _estimatedChannelMatrices;
public:
    WithChannelEstimationParticleAddon(ChannelMatrixEstimator * channelMatrixEstimator, uint trajectorylength);
    WithChannelEstimationParticleAddon(std::vector <ChannelMatrixEstimator *> channelMatrixEstimators, uint trajectorylength);
    ~WithChannelEstimationParticleAddon();
    
    WithChannelEstimationParticleAddon(const WithChannelEstimationParticleAddon& withChannelEstimationParticleAddon);

    MatrixXd getChannelMatrix_eigen(int iChannelOrder,int n) const
    {
        return _estimatedChannelMatrices[iChannelOrder][n];
    }

    void setChannelMatrix(int iChannelOrder,int n,const MatrixXd &matrix)
    {
            _estimatedChannelMatrices[iChannelOrder][n] = matrix;
    }

    ChannelMatrixEstimator *getChannelMatrixEstimator(int iChannelOrder) const
    { 
        return _channelMatrixEstimators[iChannelOrder];
    }

    int nChannelMatrixEstimators() const {return _channelMatrixEstimators.size();}

};

#endif
