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
#ifndef OCTAVE_H
#define OCTAVE_H

/**
    @author Manu <manu@rustneversleeps>
*/

#include <iostream>
#include <fstream>
#include <stdint.h>
#include <iomanip>
#include <vector>
#include <types.h>
#include <exceptions.h>

class Octave{
  
public:
    static void eigenToOctaveFileStream(const MatrixXd &A,string name,std::ofstream &f);
    static void eigenToOctaveFileStream(const std::vector<std::vector<MatrixXd> > &matrices,string name,std::ofstream &f);
    static void eigenToOctaveFileStream(const std::vector<std::vector<std::vector<MatrixXd> > > &matrices,string name,std::ofstream &f);
	static void eigenToOctaveFileStream(const std::vector<std::vector<std::vector<std::vector<MatrixXd> > > > &matrices,string name,std::ofstream &f);
    static void eigenToOctaveFileStream(const std::vector<MatrixXd> &matrices,string name,std::ofstream &f);
    template<class T> static void toOctaveFileStream(T scalar,string name,std::ofstream &f);
    static void stringsVectorToOctaveFileStream(const std::vector<string> &strings,string name,std::ofstream &f);
    template<class T> static void toOctaveFileStream(const std::vector<T> &vector,string name,std::ofstream &f);
	template<class T> static void toOctaveFileStream(const std::vector<std::vector <T> > &matrix,string name,std::ofstream &f);
	template<class T> static void toOctaveFileStream(const std::vector<std::vector<std::vector <T> > >&matrix,string name,std::ofstream &f);
};

#endif
