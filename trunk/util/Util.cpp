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

// #define DEBUG

using namespace std;

void Util::Add(const tMatrix& A,const tMatrix& B,tMatrix& C,double alpha,double beta)
{
	#ifdef DEBUG
		cout << "A.rows(): " << A.rows() << " A.cols(): " << A.cols() << " B.rows(): " << B.rows() << " B.cols(): " << B.cols() << endl;
	#endif

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

void Util::Mult(const tVector &a,const tVector &b,tMatrix &C,double alpha)
{
//     if(a.size()!=b.size() || a.size()!=C.rows() || C.rows()!=C.cols())
	if(C.rows()!=a.size() || C.cols()!=b.size())
        throw RuntimeException("Util::Mult: Resultant matrix dimensions are wrong.");

    int j;
    for(int i=0;i<a.size();i++)
    {
        for(j=0;j<b.size();j++)
            C(i,j) = alpha*a(i)*b(j);
    }
}

void Util::Transpose(const tMatrix &A,tMatrix &B)
{
	if(A.cols()!=B.rows())
		throw RuntimeException("Util::Transpose: Matrix dimensions are wrong.");

	int j;
	for(int i=0;i<A.rows();i++)
		for(j=0;j<A.cols();j++)
			B(j,i) = A(i,j);
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
	if(vector.size()> (rows*cols))
		throw RuntimeException("Util::ToMatrix: The length of the vector is greater than rows by cols.");

	tMatrix matrix = LaGenMatDouble::zeros(rows,cols);

	if(order==rowwise)
		for(uint iVector=vector.size();iVector--;)
			matrix(iVector/cols,iVector%cols) = vector(iVector);
	else
		for(uint iVector=vector.size();iVector--;)
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

	for(uint i=v.size();i--;)
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

void Util::Min(const tVector &v,int &index)
{
	double min = v(0);
	index = 0;
	for(int i=1;i<v.size();i++)
		if(v(i)<min)
		{
			min = v(i);
			index = i;
		}
}

double Util::SquareError(const tMatrix &A,const tMatrix &B)
{
    if(A.cols()!=B.cols() || A.rows()!=B.rows())
        throw IncompatibleOperandsException("Util::SquareError: matrix dimensions are different.");

    double res = 0.0;
    int j;
    for(int i=0;i<A.rows();i++)
        for(j=0;j<A.cols();j++)
            res += (A(i,j)-B(i,j))*(A(i,j)-B(i,j));
    return res;
}

void Util::Print(const tMatrix &A)
{
    int j;
    for(int i=0;i<A.rows();i++)
    {
        for(j=0;j<A.cols();j++)
            cout << setprecision(6) << setw(12) << left << A(i,j);
        cout << endl;
    }
}

void Util::MatrixToStream(tMatrix A,string name,ofstream &f)
{
    f << "# name: "<< name << endl <<"# type: matrix" << endl << "# rows: " << A.rows() << endl << "# columns: " << A.cols() << endl;

    for(int i=0;i<A.rows();i++)
    {
        for(int j=0;j<A.cols();j++)
            f << A(i,j) << " ";
        f << endl;
    }
}

void Util::MatricesVectorToStream(vector<tMatrix> matrices,string name,ofstream &f)
{
	if(matrices.size()==0 || matrices[0].rows()==0 || matrices[0].rows()==0)
		return;

	f << "# name: "<< name << endl <<"# type: matrix" << endl << "# ndims: 3" << endl << " " << (matrices.at(0)).rows() << " " << (matrices.at(0)).cols() << " " << matrices.size() << endl;

	int i,j;
	for(uint iMatrix=0;iMatrix<matrices.size();iMatrix++)
		for(j=0;j<(matrices.at(iMatrix)).cols();j++)
			for(i=0;i<(matrices.at(iMatrix)).rows();i++)
				f << " " << (matrices.at(iMatrix))(i,j) << endl;
}

template<class T> void Util::ScalarToStream(T scalar,string name,ofstream &f)
{
    f << "# name: "<< name << endl <<"# type: scalar" << endl << scalar << endl;
}

template void Util::ScalarToStream(int scalar,string name,ofstream &f);
template void Util::ScalarToStream(double scalar,string name,ofstream &f);

void Util::StringsVectorToStream(std::vector<string> strings,string name,ofstream &f)
{
	if(strings.size()==0)
		return;

	int j;

    f << "# name: "<< name << endl <<"# type: string" << endl << "# elements: " << strings.size() << endl;

	uint iMax = 0;
	uint max = strings[0].length();
	for(uint i=0;i<strings.size();i++)
		if(strings[i].length()>max)
		{
			iMax = i;
			max = strings[i].length();
		}
	for(uint i=0;i<strings.size();i++)
	{
		f << "# length: " << max << endl;
		f << strings[i];

		// paddling with spaces
		for(j=max-strings[i].length();j>0;j--)
			f << " ";
		f << endl;
	}
}

template<class T> void Util::ScalarsVectorToStream(std::vector<T> vector,string name,ofstream &f)
{
    f << "# name: "<< name << endl <<"# type: matrix" << endl << "# rows: " << "1" << endl << "# columns: " << vector.size() << endl;

	for(uint i=0;i<vector.size();i++)
		f << vector[i] << " ";
	f << endl;
}
template void Util::ScalarsVectorToStream(std::vector<double> vector,string name,ofstream &f);
template void Util::ScalarsVectorToStream(std::vector<int> vector,string name,ofstream &f);

template<class T> T Util::Max(const std::vector<T> &vector)
{
	int iMax = 0;
	T max = vector[0];

	for(uint i=1;i<vector.size();i++)
		if(vector[i]>max)
		{
			iMax = i;
			max = vector[i];
		}
	return iMax;
}
template int Util::Max(const std::vector<int> &vector);

template<class T> T Util::Sum(const std::vector<T> &vector)
{
	T sum = vector[0];

	for(uint i=1;i<vector.size();i++)
		sum += vector[i];
	return sum;
}
template int Util::Sum(const std::vector<int> &vector);

void Util::ElementByElementDiv(const tMatrix &A,const tMatrix &B,tMatrix &C)
{
	if(A.rows()!=B.rows() || A.cols()!=B.cols())
		throw RuntimeException("Util::ElementByElementDiv: Matrices can't be divided element by element.");

	int j;
	for(int i=0;i<A.rows();i++)
		for(j=0;j<A.cols();j++)
			C(i,j) = A(i,j)/B(i,j);
}

template<class T> void Util::Print(const std::vector<T> &vector)
{
	cout << "[";
	for(uint i=0;i<vector.size()-1;i++)
		cout << vector[i] << ",";
	cout << vector[vector.size()-1] << "]" << endl;
}
template void Util::Print(const std::vector<int> &vector);
template void Util::Print(const std::vector<double> &vector);

void Util::ShiftUp(tVector &v,int n)
{
	if(n>=v.size())
		throw RuntimeException("Util::ShiftUp: vector is too short for this shift.");

	for(int i=0;i<v.size()-n;i++)
		v(i) = v(i+n);
}
