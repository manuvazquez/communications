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
#include "ParticleWithChannelEstimation.h"

using namespace std;

ParticleWithChannelEstimation::ParticleWithChannelEstimation(double weight, int symbolVectorLength, int nTimeInstants
,ChannelMatrixEstimator *channelMatrixEstimator): Particle(weight, symbolVectorLength, nTimeInstants)
// ,_channelMatrixEstimator(channelMatrixEstimator),_estimatedChannelMatrices(new tMatrix[_nTimeInstants])
{
    _nChannelMatrixEstimators = 1;
    _channelMatrixEstimators = new ChannelMatrixEstimator*[_nChannelMatrixEstimators];
    _channelMatrixEstimators[0] = channelMatrixEstimator;

    _estimatedChannelMatrices = new tMatrix*[_nChannelMatrixEstimators];
    _estimatedChannelMatrices[0] = new tMatrix[_nTimeInstants];
}

// ParticleWithChannelEstimation::ParticleWithChannelEstimation(double weight, int symbolVectorLength, int nTimeInstants,vector <ChannelMatrixEstimator *> channelMatrixEstimators)
// {
//     _nChannelMatrixEstimators = channelMatrixEstimators.size();
//
// }

ParticleWithChannelEstimation::ParticleWithChannelEstimation(const ParticleWithChannelEstimation &particle):Particle(particle)
// ,_channelMatrixEstimator((particle._channelMatrixEstimator)->Clone()),_estimatedChannelMatrices(new tMatrix[_nTimeInstants])
,_nChannelMatrixEstimators(particle._nChannelMatrixEstimators)
{
    _channelMatrixEstimators = new ChannelMatrixEstimator*[particle._nChannelMatrixEstimators];
    _estimatedChannelMatrices = new tMatrix*[particle._nChannelMatrixEstimators];
    for(int iChannelMatrixEstimator=0;iChannelMatrixEstimator<particle._nChannelMatrixEstimators;iChannelMatrixEstimator++)
    {
        _channelMatrixEstimators[iChannelMatrixEstimator] = particle._channelMatrixEstimators[iChannelMatrixEstimator]->Clone();

        _estimatedChannelMatrices[iChannelMatrixEstimator] = new tMatrix[_nTimeInstants];
        for(int i=0;i<_nTimeInstants;i++)
            _estimatedChannelMatrices[iChannelMatrixEstimator][i] = particle._estimatedChannelMatrices[iChannelMatrixEstimator][i];
    }
// 	for(int i=0;i<_nTimeInstants;i++)
// 		_estimatedChannelMatrices[i] = (particle._estimatedChannelMatrices)[i];
}

ParticleWithChannelEstimation::~ParticleWithChannelEstimation()
{
// 	delete _channelMatrixEstimator;
// 	delete[] _estimatedChannelMatrices;

    for(int i=0;i<_nChannelMatrixEstimators;i++)
    {
        delete _channelMatrixEstimators[i];
        delete[] _estimatedChannelMatrices[i];
    }
    delete[] _channelMatrixEstimators;
    delete[] _estimatedChannelMatrices;
}

// void ParticleWithChannelEstimation::operator=(const ParticleWithChannelEstimation &particle)
// {
// 	if(this!=&particle)
// 	{
// 		Particle::operator =(particle);
//
// 		delete _channelMatrixEstimator;
// 		_channelMatrixEstimator = (particle._channelMatrixEstimator)->Clone();
//
// 		delete[] _estimatedChannelMatrices;
// 		_estimatedChannelMatrices = new tMatrix[_nTimeInstants];
// 		for(int i=0;i<_nTimeInstants;i++)
// 			_estimatedChannelMatrices[i] = (particle._estimatedChannelMatrices)[i];
// 	}
// }

ParticleWithChannelEstimation *ParticleWithChannelEstimation::Clone()
{
	return new ParticleWithChannelEstimation(*this);
}
