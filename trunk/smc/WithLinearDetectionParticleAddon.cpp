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
#include "WithLinearDetectionParticleAddon.h"

WithLinearDetectionParticleAddon::WithLinearDetectionParticleAddon(LinearDetector *linearDetector):_linearDetectors(1)
{
    _linearDetectors[0] = linearDetector;
}

WithLinearDetectionParticleAddon::WithLinearDetectionParticleAddon(std::vector<LinearDetector *> linearDetectors):_linearDetectors(linearDetectors)
{
}

WithLinearDetectionParticleAddon::WithLinearDetectionParticleAddon(const WithLinearDetectionParticleAddon& withLinearDetectionParticleAddon):_linearDetectors(withLinearDetectionParticleAddon._linearDetectors.size())
{
    for(uint iLinearDetector=0;iLinearDetector<_linearDetectors.size();iLinearDetector++)
        _linearDetectors[iLinearDetector] = withLinearDetectionParticleAddon._linearDetectors[iLinearDetector]->Clone();
}

WithLinearDetectionParticleAddon::~WithLinearDetectionParticleAddon()
{
    for(uint iLinearDetector=0;iLinearDetector<_linearDetectors.size();iLinearDetector++)
        delete _linearDetectors[iLinearDetector];
}


