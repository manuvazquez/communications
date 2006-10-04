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
    int _L,_N,_K;
public:
    Algorithm(string name, Alphabet  alphabet,int L,int N, int K);
	virtual ~Algorithm() {};
	virtual void Run(tMatrix observations,vector<double> noiseVariances) = 0;   
    virtual void Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence) = 0;

    string GetName() {return _name;}

    /**
    * It also returns the symbol vectors corresponding to the training sequence (if exists)
    * @return a matrix whose columns are the symbol vectors detected. It might be zero (an algorithm that knows the transmitted symbols).
    */
    virtual tMatrix GetDetectedSymbolVectors() = 0;
    /**
     * 
     * @return a vector of matrices with the channel matrices estimated. The vector length might be zero (a known channel algorithm).
     */
    virtual vector<tMatrix> GetEstimatedChannelMatrices() = 0;
    double SER(const tMatrix &symbols);
    double MSE(const vector<tMatrix> &channelMatrices);

	tMatrix HsToStackedH(vector<tMatrix> matrices,int m,int d);
	tMatrix HsToStackedH(vector<tMatrix> matrices,int m)
	{
		return HsToStackedH(matrices,m,matrices.size()-1);
	}
};

#endif
