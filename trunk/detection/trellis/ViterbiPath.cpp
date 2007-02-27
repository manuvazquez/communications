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
#include "ViterbiPath.h"

ViterbiPath::ViterbiPath():_cost(0.0),_detectedSequence(NULL)
{
}

ViterbiPath::ViterbiPath(double cost,tMatrix initialSequence):_cost(cost),_detectedSequence(new tMatrix(initialSequence))
{
}

ViterbiPath::ViterbiPath(const ViterbiPath &path):_cost(path._cost),_detectedSequence(new tMatrix(*(path._detectedSequence)))
{
}

// ViterbiPath::ViterbiPath(const ViterbiPath &path, tVector newSymbolVector, double newCost)
// {
// 	_detectedSequence = new tMatrix(path._detectedSequence->rows(),path._detectedSequence->cols()+1);
//
// 	// already detected symbols are copied into the built path
// 	(*_detectedSequence)(tRange(0,path._detectedSequence->rows()-1),tRange(0,path._detectedSequence->cols()-1)).inject(*(path._detectedSequence));
//
// 	// and so the new one
// 	_detectedSequence->col(path._detectedSequence->cols()).inject(newSymbolVector);
//
// 	// the cost is updated
// 	_cost = newCost;
// }

ViterbiPath::~ViterbiPath()
{
	delete _detectedSequence;
}

void ViterbiPath::Update(const ViterbiPath &path, tVector newSymbolVector, double newCost)
{
	if(path._detectedSequence == _detectedSequence)
		throw RuntimeException("ViterbiPath::Update=: both pointers are the same. This was not meant to happen.");

	// the below code is safe, though

	tMatrix *aux = _detectedSequence;
	_detectedSequence = new tMatrix(path._detectedSequence->rows(),path._detectedSequence->cols()+1);

	// already detected symbols are copied into the new reserved matrix
	(*_detectedSequence)(tRange(0,path._detectedSequence->rows()-1),tRange(0,path._detectedSequence->cols()-1)).inject(*(path._detectedSequence));

	// and so the new one
	_detectedSequence->col(path._detectedSequence->cols()).inject(newSymbolVector);

	// the cost is updated
	_cost = newCost;

	delete aux;
}

void ViterbiPath::Ref(const ViterbiPath &path)
{
	if(path._detectedSequence!=_detectedSequence)
	{
		delete _detectedSequence;
		_detectedSequence = path._detectedSequence;
	}
	_cost = path._cost;
}

void ViterbiPath::operator=(const ViterbiPath &path)
{
	if(path._detectedSequence == _detectedSequence)
		throw RuntimeException("ViterbiPath::operator=: both pointers are the same. This was not meant to happen.");

	// this code is safe, though
	tMatrix *aux = _detectedSequence;
	_detectedSequence = new tMatrix(*path._detectedSequence);
	delete aux;

	_cost = path._cost;
}

void ViterbiPath::Print() const
{
	cout << "Sequence:" << endl << *(_detectedSequence) << endl;
	cout << "Cost: " << _cost << endl;
}
