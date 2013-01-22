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
#include "Octave.h"

#include <bashcolors.h>

#include <sstream>

void Octave::eigenToOctaveFileStream(const MatrixXd &A,std::string name,std::ofstream &f)
{
    f << "# name: "<< name << std::endl <<"# type: matrix" << std::endl << "# rows: " << A.rows() << std::endl << "# columns: " << A.cols() << std::endl;

    for(uint i=0;i<A.rows();i++)
    {
        for(uint j=0;j<A.cols();j++)
            f << A(i,j) << " ";
        f << std::endl;
    }
}

void Octave::eigenToOctaveFileStream(const std::vector<std::vector<MatrixXd> > &matrices,std::string name,std::ofstream &f)
{
    if(matrices.size()==0 || matrices[0].size()==0 || matrices[0][0].rows()==0 || matrices[0][0].cols()==0)
    {
        std::cout << "Octave::matricesVectorsVectorToOctaveFileStream: " << COLOR_PINK << "matrix " << name << " would be an empty matrix." << COLOR_NORMAL << std::endl;
        return;
    }

    f << "# name: "<< name << std::endl <<"# type: matrix" << std::endl << "# ndims: 4" << std::endl << " " << matrices[0][0].rows() << " " << matrices[0][0].cols() << " " << matrices[0].size() << " " << matrices.size() << std::endl;

    uint i,j,k;
    for(uint l=0;l<matrices.size();l++)
        for(k=0;k<matrices[l].size();k++)
            for(j=0;j<matrices[l][k].cols();j++)
                for(i=0;i<matrices[l][k].rows();i++)
                    f << " " << matrices[l][k](i,j) << std::endl;
}

void Octave::eigenToOctaveFileStream(const std::vector<std::vector<std::vector<MatrixXd> > > &matrices,std::string name,std::ofstream &f)
{
    if(matrices.size()==0 || matrices[0].size()==0 || matrices[0][0].size()==0 || matrices[0][0][0].rows()==0 || matrices[0][0][0].cols()==0)
    {
        std::cout << "Octave::matricesVectorsVectorsVectorToOctaveFileStream: " << COLOR_PINK << "matrix " << name << " would be an empty matrix." << COLOR_NORMAL << std::endl;
        return;
    }

    f << "# name: "<< name << std::endl <<"# type: matrix" << std::endl << "# ndims: 5" << std::endl << " " << matrices[0][0][0].rows() << " " << matrices[0][0][0].cols() << " " << matrices[0][0].size() << " " << matrices[0].size() << " " << matrices.size() << std::endl;

    uint i,j,k,l;
    for(uint m=0;m<matrices.size();m++)
        for(l=0;l<matrices[m].size();l++)
            for(k=0;k<matrices[m][l].size();k++)
                for(j=0;j<matrices[m][l][k].cols();j++)
                    for(i=0;i<matrices[m][l][k].rows();i++)
                        f << " " << matrices[m][l][k](i,j) << std::endl;
}

void Octave::eigenToOctaveFileStream(const std::vector<std::vector<std::vector<std::vector<MatrixXd> > > > &matrices,std::string name,std::ofstream &f)
{
    if(matrices.size()==0 || matrices[0].size()==0 || matrices[0][0].size()==0 || matrices[0][0][0].size()==0 || matrices[0][0][0][0].rows()==0 || matrices[0][0][0][0].cols()==0)
    {
        std::cout << "Octave::matricesVectorsVectorsVectoresVectorToOctaveFileStream: " << COLOR_PINK << "matrix " << COLOR_NORMAL << name << COLOR_PINK << " would be empty: nothing is written." << COLOR_NORMAL << std::endl;
        return;
    }

    f << "# name: "<< name << std::endl <<"# type: matrix" << std::endl << "# ndims: 6" << std::endl << " " << matrices[0][0][0][0].rows() << " " << matrices[0][0][0][0].cols() << " " << matrices[0][0][0].size() << " " << matrices[0][0].size() << " " << matrices[0].size() << " " << matrices.size() << std::endl;

    uint i,j,k,l,m;
    for(uint n=0;n<matrices.size();n++)
        for(m=0;m<matrices[n].size();m++)
            for(l=0;l<matrices[n][m].size();l++)
				for(k=0;k<matrices[n][m][l].size();k++)
					for(j=0;j<matrices[n][m][l][k].cols();j++)
						for(i=0;i<matrices[n][m][l][k].rows();i++)
							f << " " << matrices[n][m][l][k](i,j) << std::endl;
}

void Octave::eigenToOctaveFileStream(const std::vector<MatrixXd> &matrices,std::string name,std::ofstream &f)
{
    uint nRows = matrices[0].rows();
	uint nCols = matrices[0].cols();
	
    if(matrices.size()==0 || nRows==0 || nCols==0)
    {
        std::cout << "Octave::matricesVectorToOctaveFileStream: " << COLOR_PINK << "matrix " << COLOR_NORMAL << name << COLOR_PINK << " would be empty: nothing is written." << COLOR_NORMAL << std::endl;
        return;
    }

    f << "# name: "<< name << std::endl <<"# type: matrix" << std::endl << "# ndims: 3" << std::endl << " " << (matrices.at(0)).rows() << " " << (matrices.at(0)).cols() << " " << matrices.size() << std::endl;

    uint i,j,iMatrix;
    for(iMatrix=0;iMatrix<matrices.size();iMatrix++)
	{
		if(matrices[iMatrix].rows()!=nRows || matrices[iMatrix].cols()!=nCols)
			throw RuntimeException("Octave::eigenToOctaveFileStream: matrices have different sizes.");
        for(j=0;j<(matrices.at(iMatrix)).cols();j++)
            for(i=0;i<(matrices.at(iMatrix)).rows();i++)
                f << " " << (matrices.at(iMatrix))(i,j) << std::endl;
	}
}

void Octave::stringsVectorToOctaveFileStream(const std::vector<std::string> &strings,std::string name,std::ofstream &f)
{
    if(strings.size()==0)
	{
		std::cout << "Octave::stringsVectorToOctaveFileStream: " << COLOR_PINK << "vector " << COLOR_NORMAL << name << COLOR_PINK << " would be empty: nothing is written." << COLOR_NORMAL << std::endl;
        return;
	}

    int j;

    f << "# name: "<< name << std::endl <<"# type: string" << std::endl << "# elements: " << strings.size() << std::endl;

	// the maximum length of all the strings is found out
    uint max = strings[0].length();
    for(uint i=0;i<strings.size();i++)
        if(strings[i].length()>max)
            max = strings[i].length();

    for(uint i=0;i<strings.size();i++)
    {
        f << "# length: " << max << std::endl;
        f << strings[i];

        // padding with spaces
        for(j=max-strings[i].length();j>0;j--)
            f << " ";
        f << std::endl;
    }
}

template<class T> void Octave::toOctaveFileStream(T scalar,std::string name,std::ofstream &f)
{
    f << "# name: "<< name << std::endl <<"# type: scalar" << std::endl << scalar << std::endl;
}

template void Octave::toOctaveFileStream(int scalar,std::string name,std::ofstream &f);
template void Octave::toOctaveFileStream(double scalar,std::string name,std::ofstream &f);
template void Octave::toOctaveFileStream(uint scalar,std::string name,std::ofstream &f);

template<class T> void Octave::toOctaveFileStream(const std::vector<T> &vector,std::string name,std::ofstream &f)
{
    f << "# name: "<< name << std::endl <<"# type: matrix" << std::endl << "# rows: " << "1" << std::endl << "# columns: " << vector.size() << std::endl;

    for(uint i=0;i<vector.size();i++)
        f << vector[i] << " ";
    f << std::endl;
}
template void Octave::toOctaveFileStream(const std::vector<double> &vector,std::string name,std::ofstream &f);
template void Octave::toOctaveFileStream(const std::vector<int> &vector,std::string name,std::ofstream &f);
template void Octave::toOctaveFileStream(const std::vector<uint32_t> &vector,std::string name,std::ofstream &f);

template<class T> void Octave::toOctaveFileStream(const std::vector<std::vector <T> > &matrix,std::string name,std::ofstream &f)
{
	uint nRows = matrix.size();
	assert(nRows!=0);
	
	f << "# name: "<< name << std::endl <<"# type: matrix" << std::endl << "# rows: " << matrix.size() << std::endl << "# columns: " << matrix[0].size() << std::endl;
	
	uint nCols = matrix[0].size();
	
	for(uint iRow=0;iRow<nRows;iRow++)
	{
		assert(matrix[iRow].size()==nCols);
		for(uint iCol=0;iCol<nCols;iCol++)
			f << matrix[iRow][iCol] << " ";
		f << std::endl;
	}
}
template void Octave::toOctaveFileStream(const std::vector<std::vector <double> > &matrix,std::string name,std::ofstream &f);
template void Octave::toOctaveFileStream(const std::vector<std::vector <bool> > &matrix,std::string name,std::ofstream &f);

template<class T> void Octave::toOctaveFileStream(const std::vector<std::vector<std::vector <T> > >&matrix,std::string name,std::ofstream &f)
{
    if(matrix.size()==0 || matrix[0].size()==0 || matrix[0][0].size()==0)
    {
        std::cout << "Octave::scalarsVectorsVectorsVectorToOctaveFileStream: " << COLOR_PINK << "matrix " << COLOR_NORMAL << name << COLOR_PINK << " would be empty: nothing is written." << COLOR_NORMAL << std::endl;
        return;
    }

    f << "# name: "<< name << std::endl <<"# type: matrix" << std::endl << "# ndims: 3" << std::endl << " " << matrix[0].size() << " " << matrix[0][0].size() << " " << matrix.size() << std::endl;

    uint i,j,iMatrix;
    for(iMatrix=0;iMatrix<matrix.size();iMatrix++)
        for(j=0;j<matrix[iMatrix][0].size();j++)
            for(i=0;i<matrix[iMatrix].size();i++)
                f << " " << matrix[iMatrix][i][j] << std::endl;
}
template void Octave::toOctaveFileStream(const std::vector<std::vector<std::vector <uint32_t> > >&matrix,std::string name,std::ofstream &f);
template void Octave::toOctaveFileStream(const std::vector<std::vector<std::vector <bool> > >&matrix,std::string name,std::ofstream &f);

MatrixXd Octave::eigenFromOctaveFileStream(std::ifstream &f)
{
	std::vector<uint> nRowsCols(2);	
	std::string read,buf;
	
	do
		getline(f,read);
	while(read.substr(0,6).compare("# type")); // ...while the first 6 characters of the line read are not equal to "# type"
	
	// number of rows and columns
	for(uint i=0;i<2;i++)
	{
		getline(f,read);
		std::stringstream ss(read);
		
		// skip # rows/columns:
		ss >> buf; ss >> buf;
		
		ss >> nRowsCols[i];
	}
	
	MatrixXd res(nRowsCols[0],nRowsCols[1]);
	
	for(uint iRow=0;iRow<nRowsCols[0];iRow++)
	{
		getline(f,read);
		std::stringstream ss(read);
		
		for(uint iCol=0;iCol<nRowsCols[1];iCol++)
			ss >> res(iRow,iCol);
	}
		
	return res;
}