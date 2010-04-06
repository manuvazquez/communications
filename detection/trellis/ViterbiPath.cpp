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

// #define DEBUG

ViterbiPath::ViterbiPath():_nTimeInstants(0),_cost(0.0),_detectedSequence(NULL)
{
}

ViterbiPath::ViterbiPath(int nTimeInstants,double cost,MatrixXd initialSequence):_nTimeInstants(nTimeInstants),_cost(cost),_detectedSequence(new MatrixXd(initialSequence))
{
}

ViterbiPath::ViterbiPath(int nTimeInstants,double cost):_nTimeInstants(nTimeInstants),_cost(cost),_detectedSequence(NULL)
{
}

ViterbiPath::ViterbiPath(const ViterbiPath &path):_nTimeInstants(path._nTimeInstants),_cost(path._cost)
{
	if(path._detectedSequence!=NULL)
		_detectedSequence = new MatrixXd(*(path._detectedSequence));
	else
		_detectedSequence = NULL;
}

ViterbiPath::~ViterbiPath()
{
	delete _detectedSequence;
}

void ViterbiPath::update(const ViterbiPath &path, VectorXd newSymbolVector, double newCost)
{
	if(path._detectedSequence == _detectedSequence)
		throw RuntimeException("ViterbiPath::Update=: both pointers are the same. This was not meant to happen.");

	// the below code is safe, though
    MatrixXd *aux = _detectedSequence;
    _detectedSequence = new MatrixXd(path._detectedSequence->rows(),path._detectedSequence->cols()+1);

	// if the accumulated sequence is not empty
	if(path._detectedSequence->cols()>0)
		// already detected symbols are copied into the new reserved matrix
        _detectedSequence->block(0,0,path._detectedSequence->rows(),path._detectedSequence->cols()) = *(path._detectedSequence);

	// and so the new one
    _detectedSequence->col(path._detectedSequence->cols()) = newSymbolVector;

	// the cost is updated
	_cost = newCost;

	// and
	_nTimeInstants = path._nTimeInstants;

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
	_nTimeInstants = path._nTimeInstants;
}

void ViterbiPath::operator=(const ViterbiPath &path)
{
	if(path._detectedSequence == _detectedSequence)
		throw RuntimeException("ViterbiPath::operator=: both pointers are the same. This was not meant to happen.");

	// this code is safe, though
    MatrixXd *aux = _detectedSequence;
    _detectedSequence = new MatrixXd(*path._detectedSequence);
	delete aux;

	_cost = path._cost;
	_nTimeInstants = path._nTimeInstants;
}

void ViterbiPath::print() const
{
	if(this->isEmpty())
		throw RuntimeException("ViterbiPath::print(): Path is empty.");

	cout << "Sequence:" << endl << *(_detectedSequence);
	cout << "Cost: " << _cost << endl;
}
