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
    static void transpose(const tMatrix &A,tMatrix &B);
    static tVector ToVector(const tMatrix &matrix,tOrder order);
    static tMatrix ToMatrix(const tVector &vector,tOrder order,int rows,int cols);
    static tMatrix ToMatrix(const vector<double> &vector,tOrder order,uint rows,uint cols);
    static tMatrix ToMatrix(const tVector &vector,tOrder order,uint rows);
    static tMatrix append(const tMatrix &A,const tMatrix &B);
    static tMatrix verticalAppend(const tMatrix &A,const tMatrix &B);
    static tVector Normalize(const tVector &v);
    static double Sum(const tVector &v);
    static void Max(const tVector &v,int &index);
    static void Min(const tVector &v,int &index);
    static double SquareErrorPaddingWithZeros(const tMatrix &A,const tMatrix &B);
    static double SquareError(const tMatrix &A,const tMatrix &B);
    static void Print(const tMatrix &A);
    static void MatrixToOctaveFileStream(tMatrix A,string name,ofstream &f);
    template<class T> static void MatricesVectorToOctaveFileStream(vector<T> matrices,string name,ofstream &f);
    static void LongIntMatricesVectorToOctaveFileStream(vector<LaGenMatLongInt> matrices,string name,ofstream &f);
    static void MatricesVectoresVectorToOctaveFileStream(vector<vector<tMatrix> > matrices,string name,ofstream &f);
    static void MatricesVectoresVectoresVectorToOctaveFileStream(vector<vector<vector<tMatrix> > > matrices,string name,ofstream &f);
    template<class T> static void ScalarToOctaveFileStream(T scalar,string name,ofstream &f);
    static void StringsVectorToOctaveFileStream(std::vector<string> strings,string name,ofstream &f);
    template<class T> static void scalarsVectorToOctaveFileStream(std::vector<T> vector,string name,ofstream &f);
    template<class T> static int Max(const std::vector<T> &vector);
    template<class T> static void Min(const std::vector<T> &vector,int &iMin);
    template<class T> static T Sum(const std::vector<T> &vector);
    static void elementWiseDiv(const tMatrix &A,const tMatrix &B,tMatrix &C);
    static void elementWiseMult(const tMatrix &A,const tMatrix &B,tMatrix &C);    
    template<class T> static void Print(const std::vector<T> &vector);
    template<class T> static void Print(const std::vector<std::vector<T> > &matrix);
    template<class T> static void Print(const T* array,int nElements);
    static void shiftUp(tVector &v,int n);
    template<class T> static vector<vector<T> > Permutations(T *array, int nElements);
    static vector<int> SolveAmbiguity(const tMatrix &H1,const tMatrix &H2,const vector<vector<uint> > &permutations,int &iBestPermutation);
    static tMatrix applyPermutation(const tMatrix &symbols,const vector<uint> &permutation,const vector<int> &signs);
    static tMatrix cholesky(const tMatrix &matrix);
    template<class T> static void NextVector(vector<T> &vector,const vector<vector<T> > &alphabets);
    template<class T> static void HowManyTimes(const vector<T> &v,vector<int> &firstOccurrence,vector<int> &times);
    static std::vector<int> nMax(int n,const tVector &v);
    static tMatrix flipLR(const tMatrix &A);
    static tMatrix sign(const tMatrix &A);   
};

#endif
