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
#ifndef UTIL_H
#define UTIL_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <math.h>
#include <vector>
#include <types.h>
#include <exceptions.h>
#include "utilExceptions.h"
#include <lapackpp/gmd.h>


using namespace std;

enum tOrder {rowwise,columnwise};

class Util{

public:

	static void Add(const tMatrix &A,const tMatrix &B,tMatrix &C,double = 1.0,double = 1.0);
	static void Add(const tVector &a,const tVector &b,tVector &c,double = 1.0,double = 1.0);
    static void Mult(const tVector &a,const tVector &b,tMatrix &C,double = 1.0);
	static void Transpose(const tMatrix &A,tMatrix &B);
	static tVector ToVector(const tMatrix &matrix,tOrder order);
	static tMatrix ToMatrix(const tVector &vector,tOrder order,int rows,int cols);
	static tMatrix ToMatrix(const tVector &vector,tOrder order,int rows);
	static tMatrix Append(const tMatrix &A,const tMatrix &B);
	static tVector Normalize(const tVector &v);
	static double Sum(const tVector &v);
	static void Max(const tVector &v,int &index);
};

#endif
