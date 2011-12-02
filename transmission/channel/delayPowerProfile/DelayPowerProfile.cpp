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
#include "DelayPowerProfile.h"

// by default the mean of the generated coefficients is zero
double DelayPowerProfile::_generatedCoefficientsMean = 0.0;

DelayPowerProfile::DelayPowerProfile(uint nOutputs,uint nInputs):_nOutputs(nOutputs),_nInputs(nInputs)
{
}

DelayPowerProfile::~DelayPowerProfile()
{
}

MatrixXd DelayPowerProfile::generateChannelMatrix(Random &random)
{
	MatrixXd res = MatrixXd::Zero(_nOutputs,_nInputs*_tapsPowers.size());

	for(uint i=0;i<res.rows();i++)
		for(uint j=0;j<res.cols();j++)
			res(i,j) = random.randn()*sqrt(_variances(i,j)) + _means(i,j);

	return res;
}

void DelayPowerProfile::print() const
{
	for(uint i=0;i<_tapsPowers.size();i++)
		std::cout << "amplitude = " << _tapsPowers[i] << std::endl;
}

void DelayPowerProfile::GenerateMatrices()
{
	// the memory of the channel is "_amplitudes.size()"
    _means = MatrixXd::Constant(_nOutputs,_nInputs*_tapsPowers.size(),_generatedCoefficientsMean);

    _variances = MatrixXd::Zero(_nOutputs,_nInputs*_tapsPowers.size());   
	for(uint i=0;i<_tapsPowers.size();i++)
		_variances.block(0,i*_nInputs,_nOutputs,_nInputs).setConstant(_tapsPowers[i]);
}
