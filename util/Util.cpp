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
#include "Util.h"

#include <bashcolors.h>

VectorXd Util::toVector(const MatrixXd &matrix,tOrder order)
{
    uint i,nElements;

    nElements = matrix.rows()*matrix.cols();
    VectorXd vector(nElements);

    if(order==rowwise)
        for(i=0;i<nElements;i++)
            vector(i) = matrix(i/matrix.cols(),i%matrix.cols());
    else
        for(i=0;i<nElements;i++)
            vector(i) = matrix(i%matrix.rows(),i/matrix.rows());
    return vector;
}

MatrixXd Util::toMatrix(const VectorXd &vector,tOrder order,uint rows,uint cols)
{
    if(vector.size()> (rows*cols))
        throw RuntimeException("Util::toMatrix: The length of the vector is greater than rows by cols.");

    MatrixXd matrix = MatrixXd::Zero(rows,cols);

    if(order==rowwise)
        for(uint iVector=vector.size();iVector--;)
            matrix(iVector/cols,iVector%cols) = vector(iVector);
    else
        for(uint iVector=vector.size();iVector--;)
            matrix(iVector%rows,iVector/rows) = vector(iVector);
    return matrix;
}

MatrixXd Util::toMatrix(const VectorXd &vector,tOrder order,uint rows)
{
    uint remainder = vector.size() % rows;
    if(remainder!=0)
        throw RuntimeException("Util::toMatrix: resultant number of columns is not integer.");
    uint cols = vector.size()/rows;
    return toMatrix(vector,order,rows,cols);
}

VectorXd Util::normalize(const VectorXd &v)
{
    int k;

    int nElements = v.size();
    double sum = 0.0;

    for(k=0;k<nElements;k++)
        sum += v(k);

    if(sum==0)
        throw AllElementsNullException("Util::normalize: A vector of zeros can't be normalized.");

    VectorXd res(nElements);
    for(k=0;k<nElements;k++)
        res(k) = v(k)/sum;
    return res;
}

void Util::normalize(std::vector<double> &v)
{
    uint k;

    uint nElements = v.size();
    double sum = 0.0;

    for(k=0;k<nElements;k++)
        sum += v[k];

    if(sum==0)
        throw AllElementsNullException("Util::normalize: A vector of zeros can't be normalized.");

    for(k=0;k<nElements;k++)
        v[k] = v[k]/sum;
}

double Util::squareErrorPaddingWithZeros(const MatrixXd &A,const MatrixXd &B)
{
    if(A.rows()!=B.rows())
	{
		cout << "A.rows() = " << A.rows() << " B.rows() = " << B.rows() << endl;
        throw RuntimeException("Util::squareErrorPaddingWithZeros: matrix have different number of rows.");
	}

    double res = 0.0;
    int i;
	int j1=-1,j2=-1; // "j1" and "j2" are initialized to -1 in case the matrices A and B have zero rows
    for(i=0;i<A.rows();i++)
        for(j1=A.cols()-1,j2=B.cols()-1;(j1>=0 && j2>=0);j1--,j2--)
            res += (A(i,j1)-B(i,j2))*(A(i,j1)-B(i,j2));

    if(j1>=0)
    {
        for(;j1>=0;j1--)
            for(i=0;i<A.rows();i++)
                res += A(i,j1)*A(i,j1);
    }
    else if(j2>=0)
    {
        for(;j2>=0;j2--)
            for(i=0;i<B.rows();i++)
                res += B(i,j2)*B(i,j2);
    }

    return res;
}

void Util::print(const MatrixXd &A)
{
    int j;
    for(int i=0;i<A.rows();i++)
    {
        for(j=0;j<A.cols();j++)
            cout << std::setprecision(6) << std::setw(12) << std::left << A(i,j);
        cout << endl;
    }
}

template<class T> T Util::sum(const std::vector<T> &vector)
{
    T sum = vector[0];

    for(uint i=1;i<vector.size();i++)
        sum += vector[i];
    return sum;
}
template int Util::sum(const std::vector<int> &vector);
template double Util::sum(const std::vector<double> &vector);

void Util::shiftUp(VectorXd &v,uint n)
{
//     if(n>=v.size())
//         throw RuntimeException("Util::shiftUp: vector is too short for this shift.");
	
	assert(n<v.size());

    for(uint i=0;i<v.size()-n;i++)
        v(i) = v(i+n);
}

template<class T> std::vector<std::vector<T> > Util::permutations(T *array, uint nElements)
{
    vector<vector<T> > res;

    uint iPermut = 0;
    do{
        res.push_back(vector<T>(nElements));
        for(uint j=0;j<nElements;j++)
            res[iPermut][j] = array[j];
        iPermut++;
    } while(std::next_permutation(array,array+nElements));

    return res;
}
template std::vector<std::vector<int> > Util::permutations(int *array, uint nElements);
template std::vector<std::vector<uint> > Util::permutations(uint *array, uint nElements);

MatrixXd Util::applyPermutationOnRows(const MatrixXd &symbols,const vector<uint> &permutation,const vector<int> &signs)
{
    uint N = symbols.rows();
    if(permutation.size()!=N || signs.size()!=N)
	{
	  std::cout << "permutation.size() = " << permutation.size() << " signs.size() = " << signs.size() << std::endl;
	  std::cout << permutation << std::endl << signs << std::endl;
	  throw RuntimeException("Util::applyPermutationOnRows: length of the received permutation is not N.");
	}

    MatrixXd res(symbols.rows(),symbols.cols());
    for(uint i=0;i<N;i++)
    {
        res.row(i) = symbols.row(permutation[i]);
        res.row(i) *= signs[permutation[i]];
    }
    return res;
}

MatrixXd Util::applyPermutationOnColumns(const MatrixXd &symbols,const vector<uint> &permutation,const vector<int> &signs)
{
    uint N = symbols.cols();
    if(permutation.size()!=N || signs.size()!=N)
	{
	  cout << "permutation.size() = " << permutation.size() << " signs.size() = " << signs.size() << " N = " << N << endl;
	  throw RuntimeException("Util::applyPermutationOnColumns: length of the received permutation is not N.");
	}

    MatrixXd res(symbols.rows(),symbols.cols());
    for(uint i=0;i<N;i++)
    {
        res.col(i) = symbols.col(permutation[i]);
        res.col(i) *= signs[permutation[i]];
    }
    return res;
}

template<class T> vector<T> Util::applyPermutation(const vector<T> &v,const vector<uint> &permutation)
{
  if(v.size()!=permutation.size())
	throw RuntimeException("Util::applyPermutation: the size of the vector to be permuted and that of the permutation doesn't match.");
  
  vector<T> res(v.size());

  for(uint i=0;i<v.size();i++)
	res[i] = v[permutation[i]];

  return res;
}
template vector<uint> Util::applyPermutation(const vector<uint> &v,const vector<uint> &permutation);
template vector<int> Util::applyPermutation(const vector<int> &v,const vector<uint> &permutation);

template<class T> void Util::nextVector(std::vector<T> &vector,const std::vector<std::vector<T> > &alphabets)
{
    if(vector.size()!=alphabets.size())
        throw RuntimeException("Util::NextVector: number of alphabets must be equal to the number of elements of the vector.");

    bool done = false;
    int iPos = vector.size()-1;
    while(iPos>=0 && !done)
    {
        uint iAlphabet = 0;
        while(vector[iPos]!=alphabets[iPos][iAlphabet] && iAlphabet<alphabets[iPos].size())
            iAlphabet++;
        if(iAlphabet==alphabets[iPos].size())
            throw RuntimeException("Util::NextVector: symbol not belonging to the corresponding alphabet found.");
        if(iAlphabet<(alphabets[iPos].size()-1))
        {
            vector[iPos] = alphabets[iPos][iAlphabet+1];
            done = true;
        }else
            vector[iPos] = alphabets[iPos][0];
        iPos--;
    }
}
template void Util::nextVector(std::vector<double> &vector,const std::vector<std::vector<double> > &alphabets);

/**
 * It finds out how many times appear each element. \ref firstOccurrence and \ref times will be deleted
 * @param v
 * @param firstOccurrence
 * @param times
 */
template<class T> void Util::howManyTimes(const vector<T> &v,vector<uint> &firstOccurrence,vector<uint> &times)
{
    firstOccurrence.clear();
    firstOccurrence.reserve(v.size());

    times.clear();
    times.reserve(v.size());

    for(uint i=0;i<v.size();i++)
    {
        uint j;
        for(j=0;j<firstOccurrence.size();j++)
            // if v[i] has already been found
            if(v[firstOccurrence[j]]==v[i])
            {
                times[j]++;
                break;
            }
        if(j==firstOccurrence.size())
        {
            firstOccurrence.push_back(i);
            times.push_back(1);
        }
    }
}
template void Util::howManyTimes(const vector<int> &v,vector<uint> &firstOccurrence,vector<uint> &times);

vector<uint> Util::nMax(uint n,const VectorXd &v)
{
    // a vector of length the minimum between the size of the vector and n is created
    vector<uint> res(n>v.size()?v.size():n);

    vector<bool> alreadySelected(v.size(),false);

    for(uint iRes=0;iRes<res.size();iRes++)
    {
        uint index = 0;
        while(alreadySelected[index])
            index++;
        double max = v(index);
        for(uint i=index+1;i<v.size();i++)
            if(!alreadySelected[i] && v(i)>max)
            {
                max = v(i);
                index = i;
            }
        res[iRes] = index;
        alreadySelected[index] = true;
    }

    return res;
}

MatrixXd Util::flipLR(const MatrixXd &A)
{
    MatrixXd res(A.rows(),A.cols());

    for(int j=0;j<A.cols();j++)
        for(int i=0;i<A.rows();i++)
            res(i,j) = A(i,A.cols()-1-j);
    return res;
}

MatrixXd Util::sign(const MatrixXd &A)
{
    MatrixXd res(A.rows(),A.cols());
    
    for(int i=0;i<A.rows();i++)
        for(int j=0;j<A.cols();j++)
            res(i,j) = (A(i,j) > 0.0)*2.0 - 1.0;
            
    return res;
}

double Util::maxCoefficientsRatio(const MatrixXd &A)
{
  double max,min;
  
  int nElements = A.rows()*A.cols();
  double currentElement;
  
	// in case, the call to this function doesn't make sense...
  if(nElements<2)
	//...it returns a non sense result (the ratio is in absolute value)
	return -1.0;
  
  max = fabs(A(0,0));
  
  // second element regardless of the size of the matrix
  min = fabs(A(1/A.cols(),1 % A.cols()));
  
  // in case we missed...
  if(max<min)
  {
	double aux = max;
	max = min;
	min = aux;
  }
  
  // we look for the max and min coefficients in absolute value
  for(int i=2;i<nElements;i++)
  {
	currentElement = fabs(A(i/A.cols(),i % A.cols()));
	if(currentElement > max)
	  max = currentElement;
	else if(currentElement < min)
	  min = currentElement;
  }
  
  return max/min;
}

std::vector<uint> Util::computeInversePermutation(const std::vector<uint> &permutation)
{
  std::vector<uint> res(permutation.size());

  uint i,j;
  
  for(i=0;i<res.size();i++)
	for(j=0;j<permutation.size();j++)
	  if(permutation[j]==i)
		res[i] = j;

  return res;
}

template<class T> std::vector<std::vector<T> > Util::block(const std::vector<std::vector<T> > &matrix, uint iStartRow, uint iStartColumn, uint nRows, uint nCols)
{
  if(iStartRow+nRows>matrix.size() || iStartColumn+nCols>matrix[0].size())
	throw RuntimeException("Util::block: not so many rows or columns.");
  
  std::vector<std::vector<T> > res(nRows,std::vector<T>(nCols));
  
  for(uint iRow=0;iRow<nRows;iRow++)
	for(uint iCol=0;iCol<nCols;iCol++)
	  res[iRow][iCol] = matrix[iRow+iStartRow][iCol+iStartColumn];
	
  return res;
}
template std::vector<std::vector<bool> > Util::block(const std::vector<std::vector<bool> > &matrix, uint iStartRow, uint iStartColumn, uint nRows, uint nCols);
template std::vector<std::vector<uint> > Util::block(const std::vector<std::vector<uint> > &matrix, uint iStartRow, uint iStartColumn, uint nRows, uint nCols);

template<class T> std::vector<T> Util::block(const std::vector<T> &vector, uint iStart, uint n)
{
	assert(iStart+n <= vector.size());
	
	std::vector<T> res(n);
	
	for(uint i=0;i<n;i++)
		res[i] = vector[iStart+i];
	
	return res;
}
template std::vector<bool> Util::block(const std::vector<bool> &vector, uint iStart, uint n);


std::vector<bool> Util::getUsersActivityFromSymbolsVector(const VectorXd &symbolsVector)
{
  std::vector<bool> res(symbolsVector.size());
  
  for(int i=0;i<symbolsVector.size();i++)
	res[i] = isUserActive(symbolsVector(i));
  
  return res;
}

template<class T> std::vector <std::vector <T > > Util::row(const std::vector <std::vector <T > > &matrix,const uint iRow)
{
	std::vector <std::vector <T > > res(1,std::vector<T>(matrix[0].size()));
	
	for(uint i=0;i<matrix[0].size();i++)
		res[0][i] = matrix[0][i];
	
	return res;
}
template std::vector <std::vector <bool> > Util::row(const std::vector <std::vector <bool> > &matrix,const uint iRow);

std::vector<MatrixXd> Util::keepCol(const std::vector<MatrixXd> &matricesVector,const uint iCol)
{
	std::vector<MatrixXd> res(matricesVector.size());
	
	for(uint i=0;i<matricesVector.size();i++)
	{
		res[i] = MatrixXd(matricesVector[i].rows(),1);
		res[i] = matricesVector[i].col(iCol);
	}
	
	return res;
}

bool Util::areColsOrthogonal(const MatrixXd &matrix)
{
	uint nRows = matrix.rows();
	uint nCols = matrix.cols();

	for (uint iOneCol=0;iOneCol<nCols;iOneCol++)
		for (uint iAnotherCol=iOneCol+1;iAnotherCol<nCols;iAnotherCol++)
		{
			int sum = 0;
			
			for (uint i=0;i<nRows;i++)
				for (uint j=0;j<nRows;j++)
					sum += matrix(i,iOneCol)*matrix(j,iAnotherCol);
			
			if (sum!=0)
				return false;
		}

	return true;
}

bool Util::areColsDifferentAndNotOpposite(const MatrixXd &matrix)
{
	uint nCols = matrix.cols();
	
	for(uint iCol=0;iCol<nCols;iCol++)
		for(uint iAnotherCol=iCol+1;iAnotherCol<nCols;iAnotherCol++)
			if( (matrix.col(iCol)==matrix.col(iAnotherCol)) || (matrix.col(iCol)==-1.0*matrix.col(iAnotherCol)))
				return false;
	return true;
}

template<class T> std::ostream& operator<<(std::ostream &out,const std::vector<T> &vector)
{
    out << "[";
    for(uint i=0;i<vector.size()-1;i++)
        out << vector[i] << ",";
    out << vector[vector.size()-1] << "]";
	
	return out;
}
template std::ostream& operator<<(std::ostream &out,const std::vector<double> &vector);
template std::ostream& operator<<(std::ostream &out,const std::vector<int> &vector);
template std::ostream& operator<<(std::ostream &out,const std::vector<uint> &vector);
template std::ostream& operator<<(std::ostream &out,const std::vector<bool> &vector);
template std::ostream& operator<<(std::ostream &out,const std::vector<MatrixXd> &vector);


template<class T> std::ostream& operator<<(std::ostream &out,const std::vector<std::vector<T> > &matrix)
{
    out << "[\n";
    for(uint i=0;i<matrix.size();i++)
    {
        for(uint j=0;j<matrix[i].size()-1;j++)
            out << matrix[i][j] << ",";
        out << matrix[i][matrix[i].size()-1] << "\n" << endl;
    }
   out << "]\n";
   
   return out;
}
template std::ostream& operator<<(std::ostream &out,const std::vector<std::vector<bool> > &matrix);
template std::ostream& operator<<(std::ostream &out,const std::vector<std::vector<uint> > &matrix);
template std::ostream& operator<<(std::ostream &out,const std::vector<std::vector<int> > &matrix);

std::vector<uint> Util::getZeroCrossings(const std::vector<MatrixXd> &matricesVector,uint iFrom, uint length)
{
	MatrixXd lastSignsMatrix = sign(matricesVector[iFrom]);

	std::vector<uint> res;
	res.push_back(iFrom);

	uint iChannelMatrix;

	for (iChannelMatrix=iFrom+1;iChannelMatrix<iFrom+length;iChannelMatrix++)
	{
		MatrixXd currentMatrixSigns = sign(matricesVector[iChannelMatrix]);
		if (currentMatrixSigns!=lastSignsMatrix)
		{
			lastSignsMatrix = currentMatrixSigns;
			res.push_back(iChannelMatrix);
		}
	}

	res.push_back(iChannelMatrix);

	return res;
}

std::vector<uint> Util::getZeroCrossings(const std::vector<MatrixXd> &matricesVector,uint iCol, uint iFrom, uint length)
{
	VectorXd lastSignsVector = sign(matricesVector[iFrom].col(iCol));

	std::vector<uint> res;
	res.push_back(iFrom);

	uint iChannelMatrix;

	for (iChannelMatrix=iFrom+1;iChannelMatrix<iFrom+length;iChannelMatrix++)
	{
		VectorXd currentVectorSigns = sign(matricesVector[iChannelMatrix].col(iCol));
		if (currentVectorSigns!=lastSignsVector)
		{
			lastSignsVector = currentVectorSigns;
			res.push_back(iChannelMatrix);
		}
	}

	res.push_back(iChannelMatrix);

	return res;
}