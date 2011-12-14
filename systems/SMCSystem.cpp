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

SMCSystem::SMCSystem()
 : BaseSystem(),ARcoefficients(1)
{
    nParticles = 1;
//     nParticles = 200;
//     nParticles = 1000;
    resamplingRatio = 0.9;

    // back and forward smoothing
    c = 0;
    e = _d;

    // AR process parameters
    ARcoefficients[0] = 0.99999;
    ARvariance=0.0001;

    // always the same resampling criterion and algorithms
    ResamplingCriterion criterioRemuestreo(resamplingRatio);

    algoritmoRemuestreo = new ResidualResamplingAlgorithm(criterioRemuestreo);

    firstSampledChannelMatrixVariance = 0.0;
}


SMCSystem::~SMCSystem()
{
  delete algoritmoRemuestreo;
}

void SMCSystem::beforeEndingFrame()
{
    BaseSystem::beforeEndingFrame();
    Octave::toOctaveFileStream(nParticles,"nParticles",_f);
    Octave::toOctaveFileStream(resamplingRatio,"resamplingRatio",_f);
    Octave::toOctaveFileStream(ARcoefficients,"ARcoefficients",_f);
    Octave::toOctaveFileStream(ARvariance,"ARvariance",_f);
    Octave::toOctaveFileStream(c,"c",_f);
    Octave::toOctaveFileStream(e,"e",_f);
}

void SMCSystem::saveFrameResults()
{
    BaseSystem::saveFrameResults();
    Octave::toOctaveFileStream(nParticles,"nParticles",_f);
    Octave::toOctaveFileStream(resamplingRatio,"resamplingRatio",_f);
    Octave::toOctaveFileStream(ARcoefficients,"ARcoefficients",_f);
    Octave::toOctaveFileStream(ARvariance,"ARvariance",_f);
    Octave::toOctaveFileStream(c,"c",_f);
    Octave::toOctaveFileStream(e,"e",_f);
}