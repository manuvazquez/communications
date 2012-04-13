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
#include "SMCSystem.h"

// #define DEBUG

#include<bashcolors.h>

SMCSystem::SMCSystem()
 : BaseSystem()
{
	xml_node<> *thisSystemParameters = get_child(_doc.first_node(),"SMCSystem");
	
	if(!thisSystemParameters)
		throw RuntimeException("SMCSystem::SMCSystem: cannot find parameters for this system.");
	
	readParameterFromXML(thisSystemParameters,"nParticles",_nParticles);
	readParameterFromXML(thisSystemParameters,"resamplingRatio",_resamplingRatio);
	readParameterFromXML(thisSystemParameters,"c",c);
	readParameterFromXML(thisSystemParameters,"firstSampledChannelMatrixVariance",_firstSampledChannelMatrixVariance);

    // always the same resampling criterion...
    ResamplingCriterion regularResamplingCriterion(_resamplingRatio);

	// ...and algorithm
    _resamplingAlgorithm = new ResidualResamplingAlgorithm(regularResamplingCriterion);
}


SMCSystem::~SMCSystem()
{
  delete _resamplingAlgorithm;
}

void SMCSystem::saveFrameResults()
{
    BaseSystem::saveFrameResults();
    Octave::toOctaveFileStream(_nParticles,"nParticles",_f);
    Octave::toOctaveFileStream(_resamplingRatio,"resamplingRatio",_f);
    Octave::toOctaveFileStream(c,"c",_f);
	Octave::toOctaveFileStream(_firstSampledChannelMatrixVariance,"firstSampledChannelMatrixVariance",_f);
}

void SMCSystem::adjustParticlesSurvivors(uint &nParticles,uint &nSurvivors,bool particlesFromSurvivors, bool survivorsFromParticles)
{
	if(particlesFromSurvivors && survivorsFromParticles)
		throw RuntimeException("SMCSystem::adjustParticlesSurvivors: both number of particles and number of survivors cannot be adjusted.");
	
	if(particlesFromSurvivors)
	{
		std::cout << COLOR_INFO << "number of particles adjusted from " << COLOR_NORMAL << nParticles;

		// the number of particles must be the number of states of the Viterbi/PSP algorithm times that of survivors
		nParticles = (uint)pow((double)_alphabet->length()+1,_N)*nSurvivors;

		std::cout << COLOR_INFO << " to " << COLOR_NORMAL << nParticles << std::endl;
	}else if(survivorsFromParticles)
	{
		cout << "Number of survivors adjusted from " << nSurvivors;
		nSurvivors = uint(ceil(double(nParticles)/pow((double)_alphabet->length()+1,double(_N))));
		cout << " to " << nSurvivors << endl;
	}
		
}

MIMOChannel *SMCSystem::createChannel()
{
	if(!_channelClassToBeInstantiated.compare(ARchannel::getXMLname()))
		return new ARchannel(_N,_L,_m,_symbols.cols(),ARprocess(_powerProfile->generateChannelMatrix(_randomGenerator),_ARcoefficients,_ARvariance));
	else
		return BaseSystem::createChannel();
}