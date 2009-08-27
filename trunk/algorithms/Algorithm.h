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
#include <LinearDetector.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas3pp.h>

class Algorithm{
protected:
    const string _name;
    const Alphabet _alphabet;
    const int _nOutputs; /// number of outputs (observations) of the system at each time instant
    const int _Nr; /// number of receiving antennas of the system
    const int _nInputs; /// number of inputs of the system at each time instant (assumed to be equal to the number of transmitting antennas/users)
    const int _iLastSymbolVectorToBeDetected;
public:
    Algorithm(string name, Alphabet  alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected);
    virtual ~Algorithm() {};

    string getName() const {return _name;}

    virtual void run(tMatrix observations,vector<double> noiseVariances)
    {
        run(Util::lapack2eigen(observations),noiseVariances);
    }
    virtual void run(MatrixXd observations,vector<double> noiseVariances) = 0;
    
    virtual void run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
    {
        run(Util::lapack2eigen(observations),noiseVariances,Util::lapack2eigen(trainingSequence));
    }
    virtual void run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence) = 0;

    virtual tMatrix getDetectedSymbolVectors()
    {
        return Util::eigen2lapack(getDetectedSymbolVectors_eigen());
    }
    /**
    * It also returns the symbol vectors corresponding to the training sequence (if it exists)
    * @return a matrix whose columns are the symbol vectors detected. It might be zero (an algorithm that knows the transmitted symbols).
    */    
    virtual MatrixXd getDetectedSymbolVectors_eigen() = 0;
    
    virtual vector<tMatrix> getEstimatedChannelMatrices()
    {
        return Util::eigen2lapack(getEstimatedChannelMatrices_eigen());
    }
    /**
     *
     * @return a vector of matrices with the channel matrices estimated. The vector length might be zero (a known channel algorithm).
     */    
    virtual vector<MatrixXd> getEstimatedChannelMatrices_eigen() = 0;


    virtual bool performsChannelOrderAPPestimation() const { return false;}

    double MSE(const vector<tMatrix> &channelMatrices)
    {
        return MSE(Util::lapack2eigen(channelMatrices));
    }
    double MSE(const vector<MatrixXd> &channelMatrices);

    tVector substractKnownSymbolsContribution(const vector<tMatrix> &matrices,int m,int c,int d,const tVector &observations,const tMatrix &symbolVectors)
    {
        return Util::eigen2lapack(substractKnownSymbolsContribution(Util::lapack2eigen(matrices),m,c,d,Util::lapack2eigen(observations),Util::lapack2eigen(symbolVectors)));
    }
    VectorXd substractKnownSymbolsContribution(const vector<MatrixXd> &matrices,int m,int c,int d,const VectorXd &observations,const MatrixXd &symbolVectors);

    tMatrix channelMatrices2stackedChannelMatrix(vector<tMatrix> matrices,int m,int start,int d)
    {
        return Util::eigen2lapack(channelMatrices2stackedChannelMatrix(Util::lapack2eigen(matrices),m,start,d));
    }
    tMatrix channelMatrices2stackedChannelMatrix(vector<tMatrix> matrices,int m)
    {
        return Util::eigen2lapack(channelMatrices2stackedChannelMatrix(Util::lapack2eigen(matrices),m,0,matrices.size()-1));
    }
    tMatrix channelMatrices2stackedChannelMatrix(vector<tMatrix> matrices,int m,int d)
    {
        return Util::eigen2lapack(channelMatrices2stackedChannelMatrix(Util::lapack2eigen(matrices),m,0,d));
    }
    
    MatrixXd channelMatrices2stackedChannelMatrix(vector<MatrixXd> matrices,int m,int start,int d);
    MatrixXd channelMatrices2stackedChannelMatrix(vector<MatrixXd> matrices,int m)
    {
        return channelMatrices2stackedChannelMatrix(matrices,m,0,matrices.size()-1);
    }
    MatrixXd channelMatrices2stackedChannelMatrix(vector<MatrixXd> matrices,int m,int d)
    {
        return channelMatrices2stackedChannelMatrix(matrices,m,0,d);
    }    
};

#endif
