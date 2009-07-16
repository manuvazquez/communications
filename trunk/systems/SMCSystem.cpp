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
    nParticles = 192;
    nParticles = 500;
//     nParticles = 1;    
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
#ifdef DEBUG
	cout << "en SMCSystem::BeforeEndingFrame" << endl;
#endif
    BaseSystem::BeforeEndingFrame(iFrame);
    Util::scalarToOctaveFileStream(nParticles,"nParticles",f);
    Util::scalarToOctaveFileStream(resamplingRatio,"resamplingRatio",f);
    Util::scalarsVectorToOctaveFileStream(ARcoefficients,"ARcoefficients",f);
    Util::scalarToOctaveFileStream(ARvariance,"ARvariance",f);
    Util::scalarToOctaveFileStream(c,"c",f);
    Util::scalarToOctaveFileStream(e,"e",f);
}
