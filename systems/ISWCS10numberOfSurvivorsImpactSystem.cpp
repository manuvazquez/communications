/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#include "ISWCS10numberOfSurvivorsImpactSystem.h"

#include <OneChannelOrderPerOutputSMCAlgorithm.h>

ISWCS10numberOfSurvivorsImpactSystem::ISWCS10numberOfSurvivorsImpactSystem():ISWCS10System()
{
	nParticlesStudied.push_back(4);
	nParticlesStudied.push_back(16);
	nParticlesStudied.push_back(128);
	nParticlesStudied.push_back(1024);
	
	
//	// 2-1-1
//	_subchannelOrders = std::vector<uint>(3,1);
//	_subchannelOrders[0] = 2;

// 	// 3-1-1
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[0] = 3;

// 	// 4-1-1
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[0] = 4;

// 	// 5-1-1
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[0] = 5;

	// 4-1-1
	_subchannelOrders = std::vector<uint>(3,1);
	_subchannelOrders[0] = 4;
	_subchannelOrders[1] = 1;
	_subchannelOrders[2] = 1;

// 	// 4-1-2
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[0] = 4;
// 	_subchannelOrders[1] = 1;
// 	_subchannelOrders[2] = 2;

// 	// 4-1-3
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[0] = 4;
// 	_subchannelOrders[1] = 1;
// 	_subchannelOrders[2] = 3;

// 	// 4-1-4
// 	_subchannelOrders = std::vector<uint>(3,1);
// 	_subchannelOrders[0] = 4;
// 	_subchannelOrders[1] = 1;
// 	_subchannelOrders[2] = 4;
}

void ISWCS10numberOfSurvivorsImpactSystem::addAlgorithms()
{
	// this is needed since the method below initializes some counters needed to deal with algorithms estimating channel orders
	ChannelOrderEstimationSystem::addAlgorithms();
	
	char algorithmName[ALGORITHM_NAME_MAX_LENGTH];

	for(uint inParticles=0;inParticles<nParticlesStudied.size();inParticles++)
	{
		sprintf(algorithmName,"MLSD-m P = %d",nParticlesStudied[inParticles]);
		_algorithms.push_back(new OneChannelOrderPerOutputSMCAlgorithm(algorithmName,*_alphabet,_L,_L,_N,_iLastSymbolVectorToBeDetected,kalmanChannelEstimators,_preamble,_preamble.cols(),_d,nParticlesStudied[inParticles],bestParticlesResamplingAlgorithm));
	}
}

void ISWCS10numberOfSurvivorsImpactSystem::saveFrameResults()
{
	ISWCS10System::saveFrameResults();
	Octave::toOctaveFileStream(nParticlesStudied,"nParticlesStudied",_f);
}