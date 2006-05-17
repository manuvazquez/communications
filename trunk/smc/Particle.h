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
#ifndef PARTICLE_H
#define PARTICLE_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <types.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>

using namespace std;

class Particle{
protected:
	double _weight;
	int _symbolVectorLength,_nTimeInstants;
	tMatrix _symbolVectors;
	tRange _rAllSymbolRows;
public:
    Particle(double weight,int symbolVectorLength,int nTimeInstants);
    virtual ~Particle();

	int TrajectoryLength() { return _nTimeInstants;}

	double GetWeight() { return _weight;}
	void SetWeight(double weight) { _weight = weight;}

	tMatrix GetAllSymbolVectors() { return _symbolVectors;}
	tVector GetSymbolVector(int n) { return _symbolVectors.col(n);}
	void SetSymbolVector(int n,const tVector &v) { _symbolVectors.col(n).inject(v);}
	void SetSymbolVector(int n,const vector<tSymbol> &v)
	{
		for(int i=0;i<_symbolVectorLength;i++)
			_symbolVectors(i,n) = v[i];
	}

	tMatrix GetSymbolVectors(const tRange &range) { return _symbolVectors(_rAllSymbolRows,range);}
	tMatrix GetSymbolVectors(int a,int b) { return _symbolVectors(_rAllSymbolRows,tRange(a,b));}
	void SetSymbolVectors(const tRange &range,const tMatrix &symbolVectors) { _symbolVectors(_rAllSymbolRows,range).inject(symbolVectors);}
	void SetSymbolVectors(int a,int b,const tMatrix &symbolVectors) { _symbolVectors(_rAllSymbolRows,tRange(a,b)).inject(symbolVectors);}

	void Print() { std::cout << _symbolVectors << std::endl << "peso = " << _weight << std::endl;}

	void operator=(const Particle &particle);
	virtual Particle *Clone();
};

#endif
