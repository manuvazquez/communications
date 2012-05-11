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

#include <iostream>
#include <stdint.h>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <types.h>
#include <exceptions.h>
#include "utilExceptions.h"

enum tOrder {rowwise,columnwise};

using std::vector;
using std::cout;
using std::endl;

class Util{
  
public:
    static VectorXd toVector(const MatrixXd &matrix,tOrder order);
    static MatrixXd toMatrix(const VectorXd &vector,tOrder order,uint rows,uint cols);
    static MatrixXd toMatrix(const VectorXd &vector,tOrder order,uint rows);
    static VectorXd normalize(const VectorXd &v);
    static void normalize(std::vector<double> &v);    
    static double squareErrorPaddingWithZeros(const MatrixXd &A,const MatrixXd &B);
    static void print(const MatrixXd &A);
    template<class T> static T sum(const std::vector<T> &vector);
    static void shiftUp(VectorXd &v,uint n);
    template<class T> static std::vector<std::vector<T> > permutations(T *array, uint nElements);
    static MatrixXd applyPermutationOnRows(const MatrixXd &symbols,const vector<uint> &permutation,const vector<int> &signs);
	static MatrixXd applyPermutationOnColumns(const MatrixXd &symbols,const vector<uint> &permutation,const vector<int> &signs);
	template<class T> static vector<T> applyPermutation(const vector<T> &v,const vector<uint> &permutation);			//
    template<class T> static void nextVector(std::vector<T> &vector,const std::vector<std::vector<T> > &alphabets);
    template<class T> static void howManyTimes(const vector<T> &v,vector<uint> &firstOccurrence,vector<uint> &times);
    
	
	/**
	 * @brief it returns the indexes of the \ref n maximum elements in \ref v
	 *
	 * @param n number of "maximum"'s
	 * @param v vector in which to look for them
	 * @return :vector< uint, std::allocator< uint > >
	 **/
	static vector<uint> nMax(uint n, const VectorXd& v);
    
	static MatrixXd flipLR(const MatrixXd &A);
    static MatrixXd sign(const MatrixXd &A);
	static std::vector<uint> computeInversePermutation(const std::vector<uint> &permutation);
	
	/*!
	  it returns the maximum ratio among coefficients of the matrix (in absolute value)
	  \param A the matrix to be analyzed
	  \return the ratio
	*/
	static double maxCoefficientsRatio(const MatrixXd &A);

	/*!
	  it returns a submatrix of the received matrix
	  \param matrix the original matrix
	  \param iStartRow the first row kept
	  \param iStartColumn the first column kept
	  \param nRows number of rows kept
	  \param nCols number of columns kept
	  \return the specified submatrix
	*/
	template<class T> static std::vector<std::vector<T> > block(const std::vector<std::vector<T> > &matrix, uint iStartRow, uint iStartColumn, uint nRows, uint nCols);
	
	/*!
	  it returns a subvector of the received vector
	  \param vector the original matrix
	  \param iStart the first element to be kept
	  \param n number of elements to be kept
	  \return the specified subvector
	*/
	template<class T> static std::vector<T> block(const std::vector<T> &vector, uint iStart, uint n);
	
	//! It returns if the the corresponding symbol implies the user is active
	/*!
	  \param symbol
	  \return true if the corresponding users is active
	*/
	static bool isUserActive(const tSymbol symbol) { return symbol!=0.0;}
	
	//! It takes a symbols vector and returns a vector indicating the active users
	/*!
	  \param symbolsVector a symbol vector
	  \return a vector of bools
	*/
	static std::vector<bool> getUsersActivityFromSymbolsVector(const VectorXd &symbolsVector);
	
	
	/**
	 * @brief it returns a row of matrix defined as a C++ vector of vectors
	 *
	 * @param matrix the matrix whose row is returned
	 * @param iRow the index of the row to be returned
	 * @return :vector< std::vector< T > >
	 **/
	template<class T> static std::vector <std::vector <T > > row(const std::vector <std::vector <T > > &matrix,const uint iRow);
	
	static std::vector<MatrixXd> keepCol(const std::vector<MatrixXd> &matricesVector,const uint iCol);
	
	/**
	 * @brief it returns the indexes where ANY of the columns (the whole column) changes sign
	 *
	 * @param matricesVector a vector of matrices
	 * @param iFrom the starting index
	 * @param length the number of matrices that must be checked starting from @ref iFrom
	 * @return a vector with the corresponding indexes. It always includes @ref iFrom as well as the last index (so that the instants define a partition of the total length)
	 **/
    static std::vector<uint> getZeroCrossings(const std::vector<MatrixXd> &matricesVector,uint iFrom, uint length);
	
	/**
	 * @brief it returns the indexes where a specific column (the whole column) changes sign
	 *
	 * @param matricesVector a vector of matrices
	 * @param iCol the column to track
	 * @param iFrom the starting index
	 * @param length the number of matrices that must be checked starting from @ref iFrom
	 * @return a vector with the corresponding indexes. It always includes @ref iFrom as well as the last index (so that the instants define a partition of the total length)
	 **/
	static std::vector<uint> getZeroCrossings(const std::vector<MatrixXd> &matricesVector,uint iCol, uint iFrom, uint length);
	
	static bool areColsOrthogonal(const MatrixXd &matrix);
	static bool areColsDifferentAndNotOpposite(const MatrixXd &matrix);
	
	static MatrixXd subMatrixFromVectorIndexes(const MatrixXd &matrix,const std::vector<uint> &rows, const std::vector<uint> &cols);
};

template<class T> std::ostream& operator<<(std::ostream &out,const std::vector<T> &vector);
template<class T> std::ostream& operator<<(std::ostream &out,const std::vector<std::vector<T> > &matrix);

#endif
