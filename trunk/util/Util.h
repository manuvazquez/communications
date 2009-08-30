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

#define EIGEN

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

    static tVector toVector(const tMatrix &matrix,tOrder order);
    static VectorXd toVector(const MatrixXd &matrix,tOrder order);
    static MatrixXd toMatrix(const VectorXd &vector,tOrder order,int rows,int cols);
    static MatrixXd toMatrix(const VectorXd &vector,tOrder order,uint rows);
    static VectorXd normalize(const VectorXd &v);
    static void normalize(std::vector<double> &v);    
    static double sum(const VectorXd &v);
    static double squareErrorPaddingWithZeros(const MatrixXd &A,const MatrixXd &B);
    static void print(const MatrixXd &A);
    static void matrixToOctaveFileStream(const MatrixXd A,string name,ofstream &f);
    template<class T> static void matricesVectorToOctaveFileStream(vector<T> matrices,string name,ofstream &f);
    static void matricesVectorsVectorToOctaveFileStream(vector<vector<MatrixXd> > matrices,string name,ofstream &f);
    static void matricesVectorsVectorsVectorToOctaveFileStream(vector<vector<vector<MatrixXd> > > matrices,string name,ofstream &f);
    template<class T> static void scalarToOctaveFileStream(T scalar,string name,ofstream &f);
    static void stringsVectorToOctaveFileStream(std::vector<string> strings,string name,ofstream &f);
    template<class T> static void scalarsVectorToOctaveFileStream(std::vector<T> vector,string name,ofstream &f);
    template<class T> static int max(const std::vector<T> &vector);
    template<class T> static void min(const std::vector<T> &vector,int &iMin);
    template<class T> static T sum(const std::vector<T> &vector);
    template<class T> static void print(const std::vector<T> &vector);
    template<class T> static void print(const std::vector<std::vector<T> > &matrix);
    template<class T> static void print(const T* array,int nElements);
    static void shiftUp(VectorXd &v,int n);
    template<class T> static vector<vector<T> > Permutations(T *array, int nElements);
    static MatrixXd applyPermutation(const MatrixXd &symbols,const vector<uint> &permutation,const vector<int> &signs);
    template<class T> static void nextVector(vector<T> &vector,const vector<vector<T> > &alphabets);
    template<class T> static void howManyTimes(const vector<T> &v,vector<int> &firstOccurrence,vector<int> &times);
    static std::vector<int> nMax(int n,const VectorXd &v);
    static MatrixXd flipLR(const MatrixXd &A);
    static MatrixXd sign(const MatrixXd &A);
    static MatrixXd lapack2eigen(const tMatrix &A);
    static VectorXd lapack2eigen(const tVector &v);
    static vector<MatrixXd> lapack2eigen(const vector<tMatrix> &v);
    static vector<vector<MatrixXd> > lapack2eigen(const vector<vector<tMatrix> > &v);
    static tMatrix eigen2lapack(const MatrixXd &A);
    static tVector eigen2lapack(const VectorXd &v);
    static vector<tMatrix> eigen2lapack(const vector<MatrixXd> &v);
};

#endif
