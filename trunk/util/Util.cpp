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

void Util::Add(const tMatrix& A,const tMatrix& B,tMatrix& C,double alpha,double beta)
{
	if(A.rows()!=B.rows() || A.cols()!=B.cols())
		throw RuntimeException("Matrices can't be added.");

	int i,j;
	int rows = A.rows(), cols = A.cols();
	for(i=0;i<rows;i++)
		for(j=0;j<cols;j++)
			C(i,j) = alpha*A(i,j) + beta*B(i,j);
}

void Util::Add(const tVector &a,const tVector &b,tVector &c,double alpha,double beta)
{
	int nElements = a.size();

	if(nElements!=b.size())
		throw RuntimeException("Vectors can't be added.");

	for(int i=0;i<nElements;i++)
		c(i) = alpha*a(i) + beta*b(i);
}

tVector Util::ToVector(const tMatrix &matrix,tOrder order)
{
	int i,nElements;

	nElements = matrix.rows()*matrix.cols();
	tVector vector(nElements);

	if(order==rowwise)
		for(i=0;i<nElements;i++)
			vector(i) = matrix(i/matrix.cols(),i%matrix.cols());
	else
		for(i=0;i<nElements;i++)
			vector(i) = matrix(i%matrix.rows(),i/matrix.rows());
	return vector;
}

tMatrix Util::ToMatrix(const tVector &vector,tOrder order,int rows,int cols)
{
	if(vector.size()!=rows*cols)
		throw RuntimeException("The length of the vector is not equal to rows by cols.");

	tMatrix matrix(rows,cols);

	if(order==rowwise)
		for(int iVector=vector.size();iVector--;)
			matrix(iVector/cols,iVector%cols) = vector(iVector);
	else
		for(int iVector=vector.size();iVector--;)
			matrix(iVector%rows,iVector/rows) = vector(iVector);
	return matrix;
}

tMatrix Util::ToMatrix(const tVector &vector,tOrder order,int rows)
{
	int remainder = vector.size() % rows;
	if(remainder!=0)
		throw RuntimeException("Resultant number of columns is not integer.");
	int cols = vector.size()/rows;
	return ToMatrix(vector,order,rows,cols);
}

tMatrix Util::Append(const tMatrix &A,const tMatrix &B)
{
	if(A.rows()!=B.rows())
		throw RuntimeException("Matrices have different number of rows.");

	tMatrix res(A.rows(),A.cols()+B.cols());
	int i,j;
	for(i=0;i<res.rows();i++)
	{
		for(j=0;j<A.cols();j++)
			res(i,j) = A(i,j);
		for(j=0;j<B.cols();j++)
			res(i,A.cols()+j) = B(i,j);
	}
	return res;
}

tVector Util::Normalize(const tVector &v)
{
	int k;

	int nElements = v.size();
	double sum = 0.0;
	
	for(k=0;k<nElements;k++)
		sum += v(k);

	if(sum==0)
		throw AllElementsNullException("Util::Normalize: A vector of zeros can't be normalized.");
	
	tVector res(nElements);
	for(k=0;k<nElements;k++)
		res(k) = v(k)/sum;
	return res;
}

double Util::Sum(const tVector &v)
{
	double res = 0.0;

	for(int i=v.size();i--;i>=0)
		res += v(i);
	return res;		
}

void Util::Max(const tVector &v,int &index)
{
	double max = v(0);
	index = 0;
	for(int i=1;i<v.size();i++)
		if(v(i)>max)
		{
			max = v(i);
			index = i;
		}
}
