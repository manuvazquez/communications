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
#ifndef ALGORITHM_H
#define ALGORITHM_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <string>
#include <vector>
#include <types.h>
#include <exceptions.h>
#include <Alphabet.h>
#include <MIMOChannel.h>
#include <lapackpp/gmd.h>
// #include <lapackpp/blas1pp.h>
// #include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>

class Algorithm{
protected:
	string _name;
	Alphabet _alphabet;
public:
    Algorithm(string name, Alphabet  alphabet);
	virtual ~Algorithm() {};
	virtual void Run(const tMatrix &observations,vector<double> noiseVariances) = 0;   
    virtual void Run(const tMatrix &observations,vector<double> noiseVariances, tMatrix trainingSequence) = 0;

    string GetName() {return _name;}

    /**
    * It also returns the symbol vectors corresponding to the training sequence (if exists)
    * @return 
    */
    virtual tMatrix GetDetectedSymbolVectors() = 0;
    virtual vector<tMatrix> GetEstimatedChannelMatrices() = 0;
    double SER(const tMatrix &symbols);
    double MSE(const vector<tMatrix> &channelMatrices);
};

#endif
