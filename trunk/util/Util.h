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

#include <iomanip>
#include <math.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <types.h>
#include <exceptions.h>
#include "utilExceptions.h"
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/sybmd.h>
#include <lapackpp/sybfd.h>

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
	static void Min(const tVector &v,int &index);
    static double SquareError(const tMatrix &A,const tMatrix &B);
    static void Print(const tMatrix &A);
	static void MatrixToStream(tMatrix A,string name,ofstream &f);
	static void MatricesVectorToStream(vector<tMatrix> matrices,string name,ofstream &f);
    template<class T> static void ScalarToStream(T scalar,string name,ofstream &f);
    static void StringsVectorToStream(std::vector<string> strings,string name,ofstream &f);
    template<class T> static void ScalarsVectorToStream(std::vector<T> vector,string name,ofstream &f);
    template<class T> static T Max(const std::vector<T> &vector);
    template<class T> static void Min(const std::vector<T> &vector,int &iMin);
    template<class T> static T Sum(const std::vector<T> &vector);
    static void ElementByElementDiv(const tMatrix &A,const tMatrix &B,tMatrix &C);
    template<class T> static void Print(const std::vector<T> &vector);
    static void ShiftUp(tVector &v,int n);
    template<class T> static vector<vector<T> > Permutations(T *array, int nElements);
    static vector<int> SolveAmbiguity(const tMatrix &H1,const tMatrix &H2,const vector<vector<uint> > &permutations,int &iBestPermutation);
    static tMatrix ApplyPermutation(const tMatrix &symbols,const vector<uint> &permutation,const vector<int> &signs);
    static tMatrix Cholesky(const tMatrix &matrix);
};

#endif
