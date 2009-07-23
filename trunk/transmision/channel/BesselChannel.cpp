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
#include "BesselChannel.h"
#define EPSILON 1e-10

// #define DEBUG

BesselChannel::BesselChannel(int nInputs, int nOutputs, int memory, int length, double velocity, double carrierFrequency, double T, const DelayPowerProfile &powerProfile): StillMemoryMIMOChannel(nInputs, nOutputs, memory, length)
{
	if(powerProfile.memory()!=memory)
		throw RuntimeException("BesselChannel::BesselChannel: memory is wrong.");

	int iRow,iCol,iTime;
	const double c = 3e8;

	double dopplerFrequency = velocity/(c/carrierFrequency);
	double normDopplerFrequency = T*dopplerFrequency;

	vector<double> autocorrelations(length);
	_channelMatrices = new tMatrix[length];

	for(iTime=0;iTime<length;iTime++)
		_channelMatrices[iTime] = tMatrix(nOutputs,nInputs*memory);

	vector<double> tapsVariances = powerProfile.tapsAmplitudes();
	vector<tMatrix> Ls(memory);
	tMatrix covarianceMatrix(length,length);

	for(uint iTap=0;iTap<tapsVariances.size();iTap++)
	{
		for(iTime=0;iTime<length;iTime++)
			autocorrelations[iTime] = tapsVariances[iTap]*jn(0,2.0*M_PI*normDopplerFrequency*double(iTime));

		for(iRow=0;iRow<length;iRow++)
			for(iCol=iRow+1;iCol<length;iCol++)
				covarianceMatrix(iCol,iRow) = covarianceMatrix(iRow,iCol) = autocorrelations[iCol-iRow];

		// variances are set
		for(iRow=0;iRow<length;iRow++)
			covarianceMatrix(iRow,iRow) = tapsVariances[iTap] + EPSILON;

		// cholesky decomposition
		Ls[iTap] = Util::cholesky(covarianceMatrix);
	}

	for(iRow=0;iRow<nOutputs;iRow++)
		for(iCol=0;iCol<nInputs*memory;iCol++)
		{
			tVector sample = LaGenMatDouble::zeros(length,1);

			// res = mean + L*randnMatrix(mean.size(),1,0.0,1.0)
			Blas_Mat_Vec_Mult(Ls[iCol/nInputs],StatUtil::randnMatrix(length,1,0.0,1.0),sample,1.0,1.0);

			for(iTime=0;iTime<length;iTime++)
				_channelMatrices[iTime](iRow,iCol) = sample(iTime);
		}
}



BesselChannel::~BesselChannel()
{
	delete[] _channelMatrices;
}


