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
private:
    static double _generatedCoefficientsMean;
protected:
    int _nOutputs,_nInputs;
	
	//! it keeps the power (variance) of each one of the channel taps
    std::vector<double> _tapsPowers;
    
// 	double _generatedCoefficientsMean;
    MatrixXd _means,_variances;

    void GenerateMatrices();
public:
    DelayPowerProfile(uint nOutputs,uint nInputs);

    virtual ~DelayPowerProfile();

    virtual MatrixXd generateChannelMatrix(Random &random);
    
    virtual void print() const;

    MatrixXd means() const { return _means;}
    
    MatrixXd variances() const {return _variances;}
   
    uint nInputs() { return _nInputs;}
    uint nOutputs() { return _nOutputs;}
    int memory() const { return _tapsPowers.size();}
    std::vector<double> tapsPowers() const { return _tapsPowers;}
    
    /**
	 * @brief It allows to set the mean of the coefficients that will be generated using this delay power profile
	 *
	 * @param mean new mean
	 * @return void
	 **/
	static void setCoefficientsMean(double mean) { _generatedCoefficientsMean = mean;}
	
	/**
	 * @brief It returns the mean of the coefficients that will be generated using this delay power profile
	 *
	 * @return double
	 **/
	static double getCoefficientsMean() { return _generatedCoefficientsMean;}
};

#endif
