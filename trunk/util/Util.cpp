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

// #define DEBUG13

using namespace std;

void Util::Add(const tMatrix& A,const tMatrix& B,tMatrix& C,double alpha,double beta)
{
	#ifdef DEBUG2
		cout << "A.rows(): " << A.rows() << " A.cols(): " << A.cols() << " B.rows(): " << B.rows() << " B.cols(): " << B.cols() << endl;
	#endif

	if(A.rows()!=B.rows() || A.cols()!=B.cols())
		throw RuntimeException("Util::Add: matrices can't be added.");

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
		throw RuntimeException("Util::Add: vectors can't be added.");

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

template<class T> void Util::Min(const std::vector<T> &vector,int &iMin)
{
    iMin = 0;
    T min = vector[0];

    for(uint i=1;i<vector.size();i++)
        if(vector[i]<min)
        {
            iMin = i;
            min = vector[i];
        }
}
template void Util::Min(const std::vector<double> &vector,int &iMin);

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
template void Util::Print(const std::vector<uint> &vector);
template void Util::Print(const std::vector<double> &vector);

void Util::ShiftUp(tVector &v,int n)
{
	if(n>=v.size())
		throw RuntimeException("Util::ShiftUp: vector is too short for this shift.");

	for(int i=0;i<v.size()-n;i++)
		v(i) = v(i+n);
}

template<class T> vector<vector<T> > Util::Permutations(T *array, int nElements)
{
	vector<vector<T> > res;

	int iPermut = 0;
	do{
		res.push_back(vector<T>(nElements));
		for(int j=0;j<nElements;j++)
			res[iPermut][j] = array[j];
		iPermut++;
    } while(next_permutation(array,array+nElements));

    return res;
}
template vector<vector<int> > Util::Permutations(int *array, int nElements);
template vector<vector<uint> > Util::Permutations(uint *array, int nElements);

vector<int> Util::SolveAmbiguity(const tMatrix &H1,const tMatrix &H2,const vector<vector<uint> > &permutations,int &iBestPermutation)
{
    #ifdef DEBUG13
        cout << "H1" << endl << H1 << "H2" << endl << H2;
    #endif

    if(H1.rows()!=H2.rows() || H1.cols()!=H2.cols())
    {
        cout << "H1" << endl << H1 << "H2" << endl << H2;
        throw RuntimeException("Util::SolveAmbiguity: matrices do not have the same dimensions.");
    }

    uint nColumns = H1.cols();

    if(permutations[0].size()!=nColumns)
        throw RuntimeException("Util::SolveAmbiguity: number of elements of the first permutations is not N.");

    for(uint i=0;i<permutations[0].size();i++)
        if(permutations[0][i]!=i)
            throw RuntimeException("Util::SolveAmbiguity: first permutation is not correct.");

    double errorWithoutChangingSign;
    double errorChangingSign;
    tVector errorVector(H1.rows());

    vector<vector<int> > signs(permutations.size(),vector<int>(nColumns));
    vector<double> permutationError(permutations.size(),0.0);

    for(uint iPermut=0;iPermut<permutations.size();iPermut++)
    {
        #ifdef DEBUG13
            cout << "probando permutación" << endl;
            Print(permutations[iPermut]);
        #endif

        for(uint iCol=0;iCol<permutations[iPermut].size();iCol++)
        {
            tVector col1 = H1.col(iCol);

            #ifdef DEBUG13
                cout << "columna de la 1ª matriz" << endl << col1;
            #endif

            // error without changing the sign
            tVector col2 = H2.col(permutations[iPermut][iCol]);

            #ifdef DEBUG13
                cout << "columna de la 2ª matriz" << endl << col2;
            #endif

            Add(col1,col2,errorVector,1.0,-1.0);
            errorWithoutChangingSign = Blas_Dot_Prod(errorVector,errorVector);

            // error changing the sign
            Add(col1,col2,errorVector,1.0,1.0);
            errorChangingSign = Blas_Dot_Prod(errorVector,errorVector);

            if(errorChangingSign<errorWithoutChangingSign)
            {
                signs[iPermut][iCol] = -1;
                permutationError[iPermut] += errorChangingSign;
            }else
            {
                signs[iPermut][iCol] = 1;
                permutationError[iPermut] += errorWithoutChangingSign;
            }
        }
    }

    #ifdef DEBUG2
        Print(permutationError);
    #endif

    Min(permutationError,iBestPermutation);
    return signs[iBestPermutation];
}

tMatrix Util::ApplyPermutation(const tMatrix &symbols,const vector<uint> &permutation,const vector<int> &signs)
{
	#ifdef DEBUG
		cout << "hola" << endl;
	#endif

    int N = symbols.rows();
    if(permutation.size()!=N || signs.size()!=N)
        throw RuntimeException("Util::ApplyPermutation: length of the received permutation is not N.");

    tMatrix res(symbols.rows(),symbols.cols());
    for(int i=0;i<N;i++)
    {
        res.row(i).inject(symbols.row(permutation[i]));
        res.row(i) *= signs[permutation[i]];
    }
	#ifdef DEBUG
		cout << "saliendo de Apply..." << endl;
	#endif
    return res;
}

// tMatrix Util::Cholesky(const tMatrix &matrix)
// {
// 	if(matrix.rows()!=matrix.cols())
// 		throw RuntimeException("Util::Cholesky: matrix is not square.");
//
// 	tMatrix res = LaGenMatDouble::zeros(matrix.rows(),matrix.rows());
//
// 	LaSymmBandMatDouble spdMatrix(matrix.rows(),2*matrix.rows()-1);
// 	for(int i=0;i<matrix.rows();i++)
// 		for(int j=0;j<matrix.cols();j++)
// 			spdMatrix(i,j) = matrix(i,j);
//
// 	LaSymmBandMatFactorizeIP(spdMatrix);
//
// 	for(int i=0;i<matrix.rows();i++)
// 		for(int j=0;j<matrix.cols();j++)
// 			if(j<=i)
// 				res(i,j) = spdMatrix(i,j);
//
// 	return res;
// }

tMatrix Util::Cholesky(const tMatrix &matrix)
{
  if (matrix.rows() != matrix.cols())
	throw RuntimeException("Util::Cholesky: Matrix not square");

  tMatrix L_ = LaGenMatDouble::zeros(matrix.rows(), matrix.rows());
  for (uint j = 0; j < matrix.rows(); j++)
	{
	  double d = 0;
	  for (uint k = 0; k < j; k++)
		{
		  double s = 0;
		  for (uint i = 0; i < k; i++)
			{
			  s += L_ (k, i) * L_ (j, i);
			}
		  L_ (j, k) = s = (matrix(j, k) - s) / L_ (k, k);
		  d = d + s * s;
		}
	  d = matrix(j, j) - d;
	  L_ (j, j) = sqrt (d);
	}
  return L_;
}
