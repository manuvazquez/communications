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

SMCSystem::SMCSystem()
 : BaseSystem(),ARcoefficients(1)
{
    nParticles = 30;
    resamplingRatio = 0.9;

    // back and forward smoothing
    c = 0;
    e = d;

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

void SMCSystem::BeforeEndingFrame(int iFrame)
{
    BaseSystem::BeforeEndingFrame(iFrame);
    Util::ScalarToOctaveFileStream(nParticles,"nParticles",f);
    Util::ScalarToOctaveFileStream(resamplingRatio,"resamplingRatio",f);
    Util::ScalarsVectorToOctaveFileStream(ARcoefficients,"ARcoefficients",f);
    Util::ScalarToOctaveFileStream(ARvariance,"ARvariance",f);
    Util::ScalarToOctaveFileStream(c,"c",f);
    Util::ScalarToOctaveFileStream(e,"e",f);
}