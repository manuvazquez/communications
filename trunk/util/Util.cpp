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

using namespace std;

void Util::add(const tMatrix& A,const tMatrix& B,tMatrix& C,double alpha,double beta)
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

void Util::add(const tVector &a,const tVector &b,tVector &c,double alpha,double beta)
{
    int nElements = a.size();

    if(nElements!=b.size())
        throw RuntimeException("Util::Add: vectors can't be added.");

    for(int i=0;i<nElements;i++)
        c(i) = alpha*a(i) + beta*b(i);
}

void Util::mult(const tVector &a,const tVector &b,tMatrix &C,double alpha)
{
//     if(a.size()!=b.size() || a.size()!=C.rows() || C.rows()!=C.cols())
    if(C.rows()!=a.size() || C.cols()!=b.size())
        throw RuntimeException("Util::mult: Resultant matrix dimensions are wrong.");

    int j;
    for(int i=0;i<a.size();i++)
    {
        for(j=0;j<b.size();j++)
            C(i,j) = alpha*a(i)*b(j);
    }
}

void Util::transpose(const tMatrix &A,tMatrix &B)
{
    if(A.cols()!=B.rows())
        throw RuntimeException("Util::transpose: Matrix dimensions are wrong.");

    int j;
    for(int i=0;i<A.rows();i++)
        for(j=0;j<A.cols();j++)
            B(j,i) = A(i,j);
}

tVector Util::toVector(const tMatrix &matrix,tOrder order)
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

// eigen
VectorXd Util::toVector(const MatrixXd &matrix,tOrder order)
{
    int i,nElements;

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

tMatrix Util::toMatrix(const tVector &vector,tOrder order,int rows,int cols)
{
    if(vector.size()> (rows*cols))
        throw RuntimeException("Util::toMatrix: The length of the vector is greater than rows by cols.");

    tMatrix matrix = LaGenMatDouble::zeros(rows,cols);

    if(order==rowwise)
        for(uint iVector=vector.size();iVector--;)
            matrix(iVector/cols,iVector%cols) = vector(iVector);
    else
        for(uint iVector=vector.size();iVector--;)
            matrix(iVector%rows,iVector/rows) = vector(iVector);
    return matrix;
}

// eigen
MatrixXd Util::toMatrix(const VectorXd &vector,tOrder order,int rows,int cols)
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

tMatrix Util::toMatrix(const tVector &vector,tOrder order,uint rows)
{
    int remainder = vector.size() % rows;
    if(remainder!=0)
        throw RuntimeException("Util::toMatrix: resultant number of columns is not integer.");
    int cols = vector.size()/rows;
    return toMatrix(vector,order,rows,cols);
}

MatrixXd Util::toMatrix(const VectorXd &vector,tOrder order,uint rows)
{
    int remainder = vector.size() % rows;
    if(remainder!=0)
        throw RuntimeException("Util::toMatrix: resultant number of columns is not integer.");
    int cols = vector.size()/rows;
    return toMatrix(vector,order,rows,cols);
}

tMatrix Util::toMatrix(const vector<double> &vector,tOrder order,uint rows,uint cols)
{
    if(vector.size()> (rows*cols))
        throw RuntimeException("Util::toMatrix: The length of the vector is greater than rows by cols.");

    tMatrix matrix = LaGenMatDouble::zeros(rows,cols);

    if(order==rowwise)
        for(uint iVector=vector.size();iVector--;)
            matrix(iVector/cols,iVector%cols) = vector[iVector];
    else
        for(uint iVector=vector.size();iVector--;)
            matrix(iVector%rows,iVector/rows) = vector[iVector];
    return matrix;
}

tMatrix Util::append(const tMatrix &A,const tMatrix &B)
{
    if(A.rows()!=B.rows())
        throw RuntimeException("Util::append: matrices have different number of rows.");

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

tMatrix Util::verticalAppend(const tMatrix &A,const tMatrix &B)
{
    if(A.cols()!=B.cols())
        throw RuntimeException("Util::verticalAppend: matrices have different number of cols.");

    tMatrix res(A.rows()+B.rows(),A.cols());
    int i,j;
    for(j=0;j<res.cols();j++)
    {
        for(i=0;i<A.rows();i++)
            res(i,j) = A(i,j);
        for(i=0;i<B.rows();i++)
            res(A.rows()+i,j) = B(i,j);
    }
    return res;
}

tVector Util::normalize(const tVector &v)
{
    int k;

    int nElements = v.size();
    double sum = 0.0;

    for(k=0;k<nElements;k++)
        sum += v(k);

    if(sum==0)
        throw AllElementsNullException("Util::normalize: A vector of zeros can't be normalized.");

    tVector res(nElements);
    for(k=0;k<nElements;k++)
        res(k) = v(k)/sum;
    return res;
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

double Util::sum(const tVector &v)
{
    double res = 0.0;

    for(uint i=v.size();i--;)
        res += v(i);
    return res;
}

// eigen
double Util::sum(const VectorXd &v)
{
    double res = 0.0;

    for(uint i=v.size();i--;)
        res += v(i);
    return res;
}

void Util::max(const tVector &v,int &index)
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

void Util::min(const tVector &v,int &index)
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

double Util::squareError(const tMatrix &A,const tMatrix &B)
{
    if(A.cols()!=B.cols() || A.rows()!=B.rows())
        throw IncompatibleOperandsException("Util::squareError: matrix dimensions are different.");

    double res = 0.0;
    int j;
    for(int i=0;i<A.rows();i++)
        for(j=0;j<A.cols();j++)
            res += (A(i,j)-B(i,j))*(A(i,j)-B(i,j));
    return res;
}

double Util::normalizedSquareError(const tMatrix &A,const tMatrix &B)
{
    if(A.cols()!=B.cols() || A.rows()!=B.rows())
        throw IncompatibleOperandsException("Util::normalizedSquareError: matrix dimensions are different.");

    double res = 0.0, normConst = 0.0;
    int j;
    for(int i=0;i<A.rows();i++)
        for(j=0;j<A.cols();j++)
        {
            res += (A(i,j)-B(i,j))*(A(i,j)-B(i,j));
            normConst = B(i,j)*B(i,j);
        }
    return res/normConst;
}

double Util::squareErrorPaddingWithZeros(const tMatrix &A,const tMatrix &B)
{
    if(A.rows()!=B.rows())
        throw IncompatibleOperandsException("Util::SquareError: matrix have different number of rows.");

    double res = 0.0;
    int i,j1,j2;
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

void Util::print(const tMatrix &A)
{
    int j;
    for(int i=0;i<A.rows();i++)
    {
        for(j=0;j<A.cols();j++)
            cout << setprecision(6) << setw(12) << left << A(i,j);
        cout << endl;
    }
}

void Util::matrixToOctaveFileStream(tMatrix A,string name,ofstream &f)
{
    f << "# name: "<< name << endl <<"# type: matrix" << endl << "# rows: " << A.rows() << endl << "# columns: " << A.cols() << endl;

    for(int i=0;i<A.rows();i++)
    {
        for(int j=0;j<A.cols();j++)
            f << A(i,j) << " ";
        f << endl;
    }
}

template<class T> void Util::matricesVectorToOctaveFileStream(vector<T> matrices,string name,ofstream &f)
{
    if(matrices.size()==0 || matrices[0].rows()==0 || matrices[0].cols()==0)
    {
        cout << "Matrix " << name << " would be an empty matrix." << endl;
        return;
    }

    f << "# name: "<< name << endl <<"# type: matrix" << endl << "# ndims: 3" << endl << " " << (matrices.at(0)).rows() << " " << (matrices.at(0)).cols() << " " << matrices.size() << endl;

    int i,j;
    for(uint iMatrix=0;iMatrix<matrices.size();iMatrix++)
        for(j=0;j<(matrices.at(iMatrix)).cols();j++)
            for(i=0;i<(matrices.at(iMatrix)).rows();i++)
                f << " " << (matrices.at(iMatrix))(i,j) << endl;
}
template void Util::matricesVectorToOctaveFileStream(vector<tMatrix> matrices,string name,ofstream &f);
template void Util::matricesVectorToOctaveFileStream(vector<LaGenMatLongInt> matrices,string name,ofstream &f);

void Util::matricesVectorsVectorToOctaveFileStream(vector<vector<tMatrix> > matrices,string name,ofstream &f)
{
    if(matrices.size()==0 || matrices[0].size()==0 || matrices[0][0].rows()==0 || matrices[0][0].cols()==0)
    {
        cout << "Matrix " << name << " would be an empty matrix." << endl;
        return;
    }

    f << "# name: "<< name << endl <<"# type: matrix" << endl << "# ndims: 4" << endl << " " << matrices[0][0].rows() << " " << matrices[0][0].cols() << " " << matrices[0].size() << " " << matrices.size() << endl;

    int i,j;
    uint k;
    for(uint l=0;l<matrices.size();l++)
        for(k=0;k<matrices[l].size();k++)
            for(j=0;j<matrices[l][k].cols();j++)
                for(i=0;i<matrices[l][k].rows();i++)
                    f << " " << matrices[l][k](i,j) << endl;
}

void Util::matricesVectorsVectorsVectorToOctaveFileStream(vector<vector<vector<tMatrix> > > matrices,string name,ofstream &f)
{
    if(matrices.size()==0 || matrices[0].size()==0 || matrices[0][0].size()==0 || matrices[0][0][0].rows()==0 || matrices[0][0][0].cols()==0)
    {
        cout << "Matrix " << name << " would be an empty matrix." << endl;
        return;
    }

    f << "# name: "<< name << endl <<"# type: matrix" << endl << "# ndims: 5" << endl << " " << matrices[0][0][0].rows() << " " << matrices[0][0][0].cols() << " " << matrices[0][0].size() << " " << matrices[0].size() << " " << matrices.size() << endl;

    int i,j;
    uint k,l;
    for(uint m=0;m<matrices.size();m++)
        for(l=0;l<matrices[m].size();l++)
            for(k=0;k<matrices[m][l].size();k++)
                for(j=0;j<matrices[m][l][k].cols();j++)
                    for(i=0;i<matrices[m][l][k].rows();i++)
                        f << " " << matrices[m][l][k](i,j) << endl;
}

template<class T> void Util::scalarToOctaveFileStream(T scalar,string name,ofstream &f)
{
    f << "# name: "<< name << endl <<"# type: scalar" << endl << scalar << endl;
}

template void Util::scalarToOctaveFileStream(int scalar,string name,ofstream &f);
template void Util::scalarToOctaveFileStream(double scalar,string name,ofstream &f);

void Util::stringsVectorToOctaveFileStream(std::vector<string> strings,string name,ofstream &f)
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

template<class T> void Util::scalarsVectorToOctaveFileStream(std::vector<T> vector,string name,ofstream &f)
{
    f << "# name: "<< name << endl <<"# type: matrix" << endl << "# rows: " << "1" << endl << "# columns: " << vector.size() << endl;

    for(uint i=0;i<vector.size();i++)
        f << vector[i] << " ";
    f << endl;
}
template void Util::scalarsVectorToOctaveFileStream(std::vector<double> vector,string name,ofstream &f);
template void Util::scalarsVectorToOctaveFileStream(std::vector<int> vector,string name,ofstream &f);
template void Util::scalarsVectorToOctaveFileStream(std::vector<uint32_t> vector,string name,ofstream &f);

template<class T> int Util::max(const std::vector<T> &vector)
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
template int Util::max(const std::vector<int> &vector);
template int Util::max(const std::vector<double> &vector);

template<class T> void Util::min(const std::vector<T> &vector,int &iMin)
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
template void Util::min(const std::vector<double> &vector,int &iMin);

template<class T> T Util::sum(const std::vector<T> &vector)
{
    T sum = vector[0];

    for(uint i=1;i<vector.size();i++)
        sum += vector[i];
    return sum;
}
template int Util::sum(const std::vector<int> &vector);
template double Util::sum(const std::vector<double> &vector);

void Util::elementWiseDiv(const tMatrix &A,const tMatrix &B,tMatrix &C)
{
    if(A.rows()!=B.rows() || A.cols()!=B.cols())
        throw RuntimeException("Util::elementWiseDiv: Matrices can't be divided element by element.");

    int j;
    for(int i=0;i<A.rows();i++)
        for(j=0;j<A.cols();j++)
            C(i,j) = A(i,j)/B(i,j);
}

void Util::elementWiseMult(const tMatrix &A,const tMatrix &B,tMatrix &C)
{
  if(A.rows()!=B.rows() || A.cols()!=B.cols())
    throw RuntimeException("Util::elementWiseMult: Matrices can't be multiplied element by element.");

  int j;
  for(int i=0;i<A.rows();i++)
    for(j=0;j<A.cols();j++)
      C(i,j) = A(i,j)*B(i,j);
}

template<class T> void Util::print(const std::vector<T> &vector)
{
    cout << "[";
    for(uint i=0;i<vector.size()-1;i++)
        cout << vector[i] << ",";
    cout << vector[vector.size()-1] << "]";
//     cout << vector[vector.size()-1] << "]" << endl;
}
template void Util::print(const std::vector<int> &vector);
template void Util::print(const std::vector<uint> &vector);
template void Util::print(const std::vector<double> &vector);
template void Util::print(const std::vector<bool> &vector);

template<class T> void Util::print(const std::vector<std::vector<T> > &matrix)
{
    cout << "[\n";
    for(uint i=0;i<matrix.size();i++)
    {
        for(uint j=0;j<matrix[i].size()-1;j++)
            cout << matrix[i][j] << ",";
        cout << matrix[i][matrix[i].size()-1] << "\n" << endl;
    }
   cout << "]\n";
}

template void Util::print(const std::vector<std::vector<bool> > &matrix);
template void Util::print(const std::vector<std::vector<uint> > &matrix);

template<class T> void Util::print(const T* array,int nElements)
{
    cout << "[";
    for(int i=0;i<nElements-1;i++)
        cout << array[i] << ",";
    cout << array[nElements-1] << "]" << endl;
}
template void Util::print(const int* array,int nElements);

void Util::shiftUp(tVector &v,int n)
{
    if(n>=v.size())
        throw RuntimeException("Util::ShiftUp: vector is too short for this shift.");

    for(int i=0;i<v.size()-n;i++)
        v(i) = v(i+n);
}

void Util::shiftUp(VectorXd &v,int n)
{
    if(n>=v.size())
        throw RuntimeException("Util::shiftUp: vector is too short for this shift.");

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

vector<int> Util::solveAmbiguity(const tMatrix &H1,const tMatrix &H2,const vector<vector<uint> > &permutations,int &iBestPermutation)
{
    #ifdef DEBUG13
        cout << "H1" << endl << H1 << "H2" << endl << H2;
    #endif

    if(H1.rows()!=H2.rows() || H1.cols()!=H2.cols())
    {
        cout << "H1" << endl << H1 << "H2" << endl << H2;
        throw RuntimeException("Util::solveAmbiguity: matrices do not have the same dimensions.");
    }

    uint nColumns = H1.cols();

    if(permutations[0].size()!=nColumns)
        throw RuntimeException("Util::solveAmbiguity: number of elements of the first permutations is not N.");

    for(uint i=0;i<permutations[0].size();i++)
        if(permutations[0][i]!=i)
            throw RuntimeException("Util::solveAmbiguity: first permutation is not correct.");

    double errorWithoutChangingSign;
    double errorChangingSign;
    tVector errorVector(H1.rows());

    vector<vector<int> > signs(permutations.size(),vector<int>(nColumns));
    vector<double> permutationError(permutations.size(),0.0);

    for(uint iPermut=0;iPermut<permutations.size();iPermut++)
    {
        #ifdef DEBUG13
            cout << "probando permutaci�n" << endl;
            print(permutations[iPermut]);
        #endif

        for(uint iCol=0;iCol<permutations[iPermut].size();iCol++)
        {
            tVector col1 = H1.col(iCol);

            #ifdef DEBUG13
                cout << "columna de la 1� matriz" << endl << col1;
            #endif

            // error without changing the sign
            tVector col2 = H2.col(permutations[iPermut][iCol]);

            #ifdef DEBUG13
                cout << "columna de la 2� matriz" << endl << col2;
            #endif

            add(col1,col2,errorVector,1.0,-1.0);
            errorWithoutChangingSign = Blas_Dot_Prod(errorVector,errorVector);

            // error changing the sign
            add(col1,col2,errorVector,1.0,1.0);
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
        print(permutationError);
    #endif

    min(permutationError,iBestPermutation);
    return signs[iBestPermutation];
}

tMatrix Util::applyPermutation(const tMatrix &symbols,const vector<uint> &permutation,const vector<int> &signs)
{
    #ifdef DEBUG
        cout << "hola" << endl;
    #endif

    uint N = symbols.rows();
    if(permutation.size()!=N || signs.size()!=N)
        throw RuntimeException("Util::ApplyPermutation: length of the received permutation is not N.");

    tMatrix res(symbols.rows(),symbols.cols());
    for(uint i=0;i<N;i++)
    {
        res.row(i).inject(symbols.row(permutation[i]));
        res.row(i) *= signs[permutation[i]];
    }
    #ifdef DEBUG
        cout << "saliendo de Apply..." << endl;
    #endif
    return res;
}

// tMatrix Util::cholesky(const tMatrix &matrix)
// {
//  if(matrix.rows()!=matrix.cols())
//      throw RuntimeException("Util::Cholesky: matrix is not square.");
//
//  tMatrix res = LaGenMatDouble::zeros(matrix.rows(),matrix.rows());
//
//  LaSymmBandMatDouble spdMatrix(matrix.rows(),2*matrix.rows()-1);
//  for(int i=0;i<matrix.rows();i++)
//      for(int j=i;j<matrix.cols();j++)
//          spdMatrix(j,i) = spdMatrix(i,j) = matrix(i,j);
//
//  LaSymmBandMatFactorizeIP(spdMatrix);
//
//  for(int i=0;i<matrix.rows();i++)
//      for(int j=0;j<matrix.cols();j++)
//          if(j<=i)
//              res(i,j) = spdMatrix(i,j);
//
//  return res;
// }

tMatrix Util::cholesky(const tMatrix &matrix)
{
  if (matrix.rows() != matrix.cols())
    throw RuntimeException("Util::Cholesky: Matrix not square");

  tMatrix L_ = LaGenMatDouble::zeros(matrix.rows(), matrix.rows());
  for (int j = 0; j < matrix.rows(); j++)
    {
      double d = 0;
      for (int k = 0; k < j; k++)
        {
          double s = 0;
          for (int i = 0; i < k; i++)
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

// eigen
MatrixXd Util::cholesky(const MatrixXd &matrix)
{
  if (matrix.rows() != matrix.cols())
    throw RuntimeException("Util::Cholesky: Matrix not square");

  MatrixXd L_ = MatrixXd::Zero(matrix.rows(), matrix.rows());
  for (int j = 0; j < matrix.rows(); j++)
    {
      double d = 0.0;
      for (int k = 0; k < j; k++)
        {
          double s = 0.0;
          for (int i = 0; i < k; i++)
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

template<class T> void Util::nextVector(vector<T> &vector,const vector<vector<T> > &alphabets)
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
template void Util::nextVector(vector<double> &vector,const vector<vector<double> > &alphabets);

/**
 * It finds out how many times appear each element. @param firstOccurrence and @param times will be deleted
 * @param v
 * @param firstOccurrence
 * @param times
 */
template<class T> void Util::howManyTimes(const vector<T> &v,vector<int> &firstOccurrence,vector<int> &times)
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
template void Util::howManyTimes(const vector<int> &v,vector<int> &firstOccurrence,vector<int> &times);

vector<int> Util::nMax(int n,const tVector &v)
{
    // a vector of length the minimum between the size of the vector and n is created
    vector<int> res(n>v.size()?v.size():n);

    vector<bool> alreadySelected(v.size(),false);

    for(uint iRes=0;iRes<res.size();iRes++)
    {
        int index = 0;
        while(alreadySelected[index])
            index++;
        double max = v(index);
        for(int i=index+1;i<v.size();i++)
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

vector<int> Util::nMax(int n,const VectorXd &v)
{
    // a vector of length the minimum between the size of the vector and n is created
    vector<int> res(n>v.size()?v.size():n);

    vector<bool> alreadySelected(v.size(),false);

    for(uint iRes=0;iRes<res.size();iRes++)
    {
        int index = 0;
        while(alreadySelected[index])
            index++;
        double max = v(index);
        for(int i=index+1;i<v.size();i++)
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

tMatrix Util::flipLR(const tMatrix &A)
{
    tMatrix res(A.rows(),A.cols());

    for(int j=0;j<A.cols();j++)
        for(int i=0;i<A.rows();i++)
            res(i,j) = A(i,A.cols()-1-j);
    return res;
}

// eigen
MatrixXd Util::flipLR(const MatrixXd &A)
{
    MatrixXd res(A.rows(),A.cols());

    for(int j=0;j<A.cols();j++)
        for(int i=0;i<A.rows();i++)
            res(i,j) = A(i,A.cols()-1-j);
    return res;
}

tMatrix Util::sign(const tMatrix &A)
{
    tMatrix res(A.rows(),A.cols());
    
    for(int i=0;i<A.rows();i++)
        for(int j=0;j<A.cols();j++)
            res(i,j) = (A(i,j) > 0.0)*2.0 - 1.0;
            
    return res;
}

MatrixXd Util::lapack2eigen(const tMatrix &A)
{
    int i,j;
    MatrixXd res(A.rows(),A.cols());
    
    for(i=0;i<A.rows();i++)
        for(j=0;j<A.cols();j++)
            res(i,j) = A(i,j);
            
    return res;
}
VectorXd Util::lapack2eigen(const tVector &v)
{
    int i;
    VectorXd res(v.size());
    
    for(i=0;i<v.size();i++)
        res(i) = v(i);
    
    return res;
}

vector<MatrixXd> Util::lapack2eigen(const vector<tMatrix> &v)
{
    vector<MatrixXd> res(v.size());
    
    for(uint i=0;i<v.size();i++)
        res[i] = Util::lapack2eigen(v[i]);
        
    return res;
}

vector<vector<MatrixXd> > Util::lapack2eigen(const vector<vector<tMatrix> > &v)
{
    vector<vector<MatrixXd> > res(v.size());
    
    for(uint i=0;i<v.size();i++)
        for(uint j=0;j<v[i].size();j++)
            res[i].push_back(Util::lapack2eigen(v[i][j]));
        
    return res;
}

tMatrix Util::eigen2lapack(const MatrixXd &A)
{
    int i,j;
    tMatrix res(A.rows(),A.cols());
    
    for(i=0;i<A.rows();i++)
        for(j=0;j<A.cols();j++)
            res(i,j) = A(i,j);
            
    return res;
}
tVector Util::eigen2lapack(const VectorXd &v)
{
    int i;
    tVector res(v.size());
    
    for(i=0;i<v.size();i++)
        res(i) = v(i);
    
    return res;
}

vector<tMatrix> Util::eigen2lapack(const vector<MatrixXd> &v)
{
    vector<tMatrix> res(v.size());
    
    for(uint i=0;i<v.size();i++)
        res[i] = Util::eigen2lapack(v[i]);
        
    return res;
}