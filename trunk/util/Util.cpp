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

MatrixXd Util::toMatrix(const VectorXd &vector,tOrder order,uint rows)
{
    int remainder = vector.size() % rows;
    if(remainder!=0)
        throw RuntimeException("Util::toMatrix: resultant number of columns is not integer.");
    int cols = vector.size()/rows;
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

double Util::sum(const VectorXd &v)
{
    double res = 0.0;

    for(uint i=v.size();i--;)
        res += v(i);
    return res;
}

double Util::squareErrorPaddingWithZeros(const MatrixXd &A,const MatrixXd &B)
{
    if(A.rows()!=B.rows())
        throw IncompatibleOperandsException("Util::squareError: matrix have different number of rows.");

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

void Util::print(const MatrixXd &A)
{
    int j;
    for(int i=0;i<A.rows();i++)
    {
        for(j=0;j<A.cols();j++)
            cout << setprecision(6) << setw(12) << left << A(i,j);
        cout << endl;
    }
}

void Util::matrixToOctaveFileStream(const MatrixXd A,string name,ofstream &f)
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
// template void Util::matricesVectorToOctaveFileStream(vector<LaGenMatLongInt> matrices,string name,ofstream &f);
template void Util::matricesVectorToOctaveFileStream(vector<MatrixXd> matrices,string name,ofstream &f);

void Util::matricesVectorsVectorToOctaveFileStream(vector<vector<MatrixXd> > matrices,string name,ofstream &f)
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

void Util::matricesVectorsVectorsVectorToOctaveFileStream(vector<vector<vector<MatrixXd> > > matrices,string name,ofstream &f)
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

template<class T> void Util::print(const std::vector<T> &vector)
{
    cout << "[";
    for(uint i=0;i<vector.size()-1;i++)
        cout << vector[i] << ",";
    cout << vector[vector.size()-1] << "]";
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

MatrixXd Util::applyPermutation(const MatrixXd &symbols,const vector<uint> &permutation,const vector<int> &signs)
{
    uint N = symbols.rows();
    if(permutation.size()!=N || signs.size()!=N)
        throw RuntimeException("Util::ApplyPermutation: length of the received permutation is not N.");

    MatrixXd res(symbols.rows(),symbols.cols());
    for(uint i=0;i<N;i++)
    {
        res.row(i) = symbols.row(permutation[i]);
        res.row(i) *= signs[permutation[i]];
    }
    return res;
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

//! first (inefficient) implementation of \ref maxCoefficientsRatio
double Util::maxCoefficientsRatio2(const MatrixXd &A)
{
  int iMax1,jMax1,iMax2,jMax2;
  double max = -1.0;
  double ratio;
  
  int i,j,i2,j2;
  
  for(i=0;i<A.rows();i++)
	for(j=0;j<A.cols();j++)
	{
	  for(i2=i;i2<A.rows();i2++)
		for(j2=j+1;j2<A.cols();j2++)
		{
		  ratio = fabs(A(i,j)/A(i2,j2));
		  
		  // if we have divided the smaller coefficient by the larger one
		  if(ratio<1.0)
			ratio = 1.0/ratio;
		  
		  if(ratio>max)
		  {
			max = ratio;
			iMax1 = i;
			jMax1 = j;
			iMax2 = i2;
			jMax2 = j2;
		  }
		}
	}
	
	return max;
}

double Util::maxCoefficientsRatio(const MatrixXd &A)
{
  double max,min;
  
  int nElements = A.rows()*A.cols();
  double currentElement;
  
  if(nElements<2)
	throw RuntimeException("Util::maxCoefficientsRatio: the channel has a single coefficient.");
  
  max = fabs(A(0,0));
  min = fabs(A(1/A.cols(),1 % A.cols()));
  
  // in case we missed...
  if(max<min)
  {
	double aux = max;
	max = min;
	min = aux;
  }
  
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