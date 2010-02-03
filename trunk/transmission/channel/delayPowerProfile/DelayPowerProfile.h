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
#ifndef DELAYPOWERPROFILE_H
#define DELAYPOWERPROFILE_H

/**
    @author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <vector>
#include <Random.h>
#include <Util.h>

class DelayPowerProfile{
protected:
    int _nOutputs,_nInputs;
    std::vector<double> _amplitudes;
    double _generatedCoefficientsMean;
    MatrixXd _means,_variances;

    void GenerateMatrices();
public:
    DelayPowerProfile(int nOutputs,int nInputs);

    virtual ~DelayPowerProfile();

    virtual MatrixXd generateChannelMatrix(Random &random);
    
    virtual void print() const;

    MatrixXd means() const { return _means;}
    
    MatrixXd variances() const {return _variances;}
    
    int nInputs() { return _nInputs;}
    int nOutputs() { return _nOutputs;}
    int memory() const { return _amplitudes.size();}
    std::vector<double> tapsAmplitudes() const { return _amplitudes;}
};

#endif
