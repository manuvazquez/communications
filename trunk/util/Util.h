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

#include <types.h>
#include <Random.h>
#include <lapackpp/gmd.h>

class Util{

private:
// 	static Random randomGenerator(234234);

public:
    Util() {};
    ~Util() {};

	static void Add(const tMatrix& A,const tMatrix& B,tMatrix& C,double = 1.0,double = 1.0);
// 	static void GaussianMatrix(tMatrix& matrix,double = 0,double = 1);
};

#endif
