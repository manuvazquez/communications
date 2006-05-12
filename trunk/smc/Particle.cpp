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
#include "Particle.h"

Particle::Particle(double weight,int symbolVectorLength,int nTimeInstants):_weight(weight),_symbolVectorLength(symbolVectorLength),_nTimeInstants(nTimeInstants),_symbolVectors(_symbolVectorLength,_nTimeInstants),_rAllSymbolRows(0,_symbolVectorLength-1)
{
}

Particle::~Particle()
{
}

void Particle::operator=(const Particle &particle)
{
// 	double _weight;
// 	int _symbolVectorLength,_nTimeInstants;
// 	tMatrix _symbolVectors;
// 	tRange _rAllSymbolRows;
	if(this!=&particle)
	{
		_weight = particle._weight;
		_symbolVectorLength = particle._symbolVectorLength;
		_nTimeInstants = particle._nTimeInstants;
		_symbolVectors = particle._symbolVectors;
		_rAllSymbolRows = particle._rAllSymbolRows;
	}
}

Particle *Particle::Clone()
{
	return new Particle(*this);
}
