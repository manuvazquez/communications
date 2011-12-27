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
	
	readParameterFromXML(thisSystemParameters,"nParticles",nParticles);
	readParameterFromXML(thisSystemParameters,"resamplingRatio",resamplingRatio);
	readParameterFromXML(thisSystemParameters,"c",c);
	
	xml_node<> *ARprocessNode = get_child(thisSystemParameters,"ARprocess");
	if(!ARprocessNode)
		throw RuntimeException("SMCSystem::SMCSystem: cannot find parameter \"ARprocess\"");
	readMultiValuedParameterFromXML(ARprocessNode,"coefficients",_ARcoefficients);
	readParameterFromXML(ARprocessNode,"variance",_ARvariance);
	
// // 	cout << "nParticles = " << nParticles << " resamplingRatio = " << resamplingRatio << " c = " << c << endl;
// // 	cout << "_ARcoefficients = " << endl << _ARcoefficients << endl << "_ARvariance = " << _ARvariance << endl;
// 	
// 	
//     nParticles = 1;
// //     nParticles = 200;
// //     nParticles = 1000;
//     resamplingRatio = 0.9;
// 
//     // back and forward smoothing
//     c = 0;
// 
//     // AR process parameters
// //     _ARcoefficients[0] = 0.99999;
// 	_ARcoefficients.push_back(0.99999);
//     _ARvariance=0.0001;
// 	
//     firstSampledChannelMatrixVariance = 0.0;

    // always the same resampling criterion and algorithms
    ResamplingCriterion criterioRemuestreo(resamplingRatio);

    algoritmoRemuestreo = new ResidualResamplingAlgorithm(criterioRemuestreo);
}


SMCSystem::~SMCSystem()
{
  delete algoritmoRemuestreo;
}

void SMCSystem::saveFrameResults()
{
    BaseSystem::saveFrameResults();
    Octave::toOctaveFileStream(nParticles,"nParticles",_f);
    Octave::toOctaveFileStream(resamplingRatio,"resamplingRatio",_f);
    Octave::toOctaveFileStream(_ARcoefficients,"ARcoefficients",_f);
    Octave::toOctaveFileStream(_ARvariance,"ARvariance",_f);
    Octave::toOctaveFileStream(c,"c",_f);
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