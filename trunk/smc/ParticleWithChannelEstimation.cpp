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

ParticleWithChannelEstimation::ParticleWithChannelEstimation(double weight, int symbolVectorLength, int nTimeInstants
// ,int channelMatrixRows,int channelMatrixColumns
,ChannelMatrixEstimator *channelMatrixEstimator): Particle(weight, symbolVectorLength, nTimeInstants)
// ,_channelMatrixRows(channelMatrixRows),_channelMatrixColumns(channelMatrixColumns)
,_channelMatrixEstimator(channelMatrixEstimator),_estimatedChannelMatrices(new tMatrix[_nTimeInstants])
{
}

ParticleWithChannelEstimation::ParticleWithChannelEstimation(const ParticleWithChannelEstimation &particle):Particle(particle)
// ,_channelMatrixRows(particle._channelMatrixRows),_channelMatrixColumns(particle._channelMatrixColumns)
,_channelMatrixEstimator((particle._channelMatrixEstimator)->Clone()),_estimatedChannelMatrices(new tMatrix[_nTimeInstants])
{
	for(int i=0;i<_nTimeInstants;i++)
		_estimatedChannelMatrices[i] = (particle._estimatedChannelMatrices)[i];
}

ParticleWithChannelEstimation::~ParticleWithChannelEstimation()
{
	delete _channelMatrixEstimator;
	delete[] _estimatedChannelMatrices;
}

void ParticleWithChannelEstimation::operator=(const ParticleWithChannelEstimation &particle)
{
	if(this!=&particle)
	{
		Particle::operator =(particle);

// 		_channelMatrixRows = particle._channelMatrixRows;
// 		_channelMatrixColumns = particle._channelMatrixColumns;

		delete _channelMatrixEstimator;
		_channelMatrixEstimator = (particle._channelMatrixEstimator)->Clone();

		delete[] _estimatedChannelMatrices;
		_estimatedChannelMatrices = new tMatrix[_nTimeInstants];
		for(int i=0;i<_nTimeInstants;i++)
			_estimatedChannelMatrices[i] = (particle._estimatedChannelMatrices)[i];		
	}
}

ParticleWithChannelEstimation *ParticleWithChannelEstimation::Clone()
{
	return new ParticleWithChannelEstimation(*this);
}
