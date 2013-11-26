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

#include <defines.h>
#include <string>
#include <vector>
#include <types.h>
#include <exceptions.h>
#include <Alphabet.h>
#include <Util.h>

class Algorithm{
protected:
    const std::string _name;
    const Alphabet _alphabet;
    const uint _nOutputs; /// number of outputs (observations) of the system at each time instant
    const uint _Nr; /// number of receiving antennas of the system
    const uint _nInputs; /// number of inputs of the system at each time instant (assumed to be equal to the number of transmitting antennas/users)
    const uint _iLastSymbolVectorToBeDetected;
public:
    Algorithm(std::string name, Alphabet  alphabet,uint L,uint Nr,uint N, uint iLastSymbolVectorToBeDetected);
    virtual ~Algorithm() {};

    std::string getName() const {return _name;}

    virtual void run(MatrixXd observations,std::vector<double> noiseVariances) = 0;
    
    virtual void run(MatrixXd observations,std::vector<double> noiseVariances, MatrixXd trainingSequence) = 0;

    /*!
    * It also returns the symbol vectors corresponding to the training sequence (if it exists)
    * \return a matrix whose columns are the detected symbol vectors. It might be zero (an algorithm that knows the transmitted symbols).
    */    
    virtual MatrixXd getDetectedSymbolVectors() = 0;
    
    /*!
     *
     * \return a vector of matrices with the channel matrices estimated. The vector length might be zero (a known channel algorithm).
     */    
    virtual std::vector<MatrixXd> getEstimatedChannelMatrices() = 0;
	
#ifdef SAVE_CHANNEL_ESTIMATES_VARIANCES
	virtual std::vector<MatrixXd> getChannelEstimatesVariances() const { throw RuntimeException("Algorithm::getChannelEstimatesVariances: not implemented for this algorithm."); }
#endif

    /*!
    * It also returns the symbol vectors corresponding to the training sequence (if it exists)
    * \return a matrix whose columns are the estimated symbol vectors. It might be zero (an algorithm that knows the transmitted symbols).
    */    
    virtual MatrixXd getEstimatedSymbolVectors() { throw RuntimeException("Algorithm::getEstimatedSymbolVectors: not implemented for this algorithm."); }


    virtual bool estimatesOneSingleChannelOrder() const { return false;}
    virtual bool performsSymbolsDetection() const { return true; }
    
    virtual bool performsChannelEstimation() const { return true; }
    virtual bool computesChannelEstimatesVariances() const { return false; }
    
    virtual bool estimatesOneChannelOrderPerOutput() const { return false;}
    
    virtual bool performsSymbolsEstimation() const { return false; }

    VectorXd substractKnownSymbolsContribution(const std::vector<MatrixXd> &matrices,uint m,uint d,const VectorXd &observations,const MatrixXd &symbolVectors);
};

#endif
