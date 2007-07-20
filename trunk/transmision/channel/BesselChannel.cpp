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

BesselChannel::BesselChannel(int nTx, int nRx, int memory, int length, double velocity, double carrierFrequency, double T, const DelayPowerProfile &powerProfile): StillMemoryMIMOChannel(nTx, nRx, memory, length)
{
	if(powerProfile.Memory()!=memory)
		throw RuntimeException("BesselChannel::BesselChannel: memory is wrong.");

	int iRow,iCol,iTime;
	const double c = 3e8;

	double dopplerFrequency = velocity/(c/carrierFrequency);
	double normDopplerFrequency = T*dopplerFrequency;

	vector<double> autocorrelations(length);
	_channelMatrices = new tMatrix[length];

	for(iTime=0;iTime<length;iTime++)
		_channelMatrices[iTime] = tMatrix(nRx,nTx*memory);

	vector<double> tapsVariances = powerProfile.TapsAmplitudes();
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
		Ls[iTap] = Util::Cholesky(covarianceMatrix);
	}

	for(iRow=0;iRow<nRx;iRow++)
		for(iCol=0;iCol<nTx*memory;iCol++)
		{
			tVector sample = LaGenMatDouble::zeros(length,1);

			// res = mean + L*RandnMatrix(mean.size(),1,0.0,1.0)
			Blas_Mat_Vec_Mult(Ls[iCol/nTx],StatUtil::RandnMatrix(length,1,0.0,1.0),sample,1.0,1.0);

			for(iTime=0;iTime<length;iTime++)
				_channelMatrices[iTime](iRow,iCol) = sample(iTime);
		}
}

// BesselChannel::BesselChannel(int nTx, int nRx, int memory, int length, double velocity, double carrierFrequency, double T, const ContinuousPowerProfile &powerProfile): StillMemoryMIMOChannel(nTx, nRx, memory, length)
// {
// 	int iRow,iCol,iTime;
// 	const double c = 3e8;
//
// 	double dopplerFrequency = velocity/(c/carrierFrequency);
// 	double normDopplerFrequency = T*dopplerFrequency;
//
// 	vector<double> autocorrelations(length);
//
// 	vector<tMatrix> profiledChannelMatrices(length);
// 	_channelMatrices = new tMatrix[length];
//
// 	vector<tMatrix> channelMatricesOut(length);
//
// 	for(iTime=0;iTime<length;iTime++)
// 	{
// 		profiledChannelMatrices[iTime] = tMatrix(nRx,nTx*powerProfile.Ndelays());
// 		_channelMatrices[iTime] = LaGenMatDouble::zeros(nRx,nTx*memory);
// 	}
//
// 	vector<double> powers = powerProfile.Powers();
// 	vector<tMatrix> Ls(powerProfile.Ndelays());
// 	tMatrix covarianceMatrix(length,length);
//
// #ifdef DEBUG
// 	cout << "powers.size() = " << endl << powers.size() << endl;
// 	cout << "Los retardos" << endl;
// 	Util::Print(powerProfile.Delays());
// 	cout << "Las potencias" << endl;
// 	Util::Print(powerProfile.Powers());
// #endif
//
// 	for(uint iPower=0;iPower<powers.size();iPower++)
// 	{
// 		for(iTime=0;iTime<length;iTime++)
// 			autocorrelations[iTime] = powers[iPower]*jn(0,2.0*M_PI*normDopplerFrequency*double(iTime));
//
// 		for(iRow=0;iRow<length;iRow++)
// 			for(iCol=iRow+1;iCol<length;iCol++)
// 				covarianceMatrix(iCol,iRow) = covarianceMatrix(iRow,iCol) = autocorrelations[iCol-iRow];
//
// 		// variances are set
// 		for(iRow=0;iRow<length;iRow++)
// 			covarianceMatrix(iRow,iRow) = powers[iPower] + EPSILON;
//
// 		// cholesky decomposition
// 		Ls[iPower] = Util::Cholesky(covarianceMatrix);
// 	}
//
// 	for(iRow=0;iRow<nRx;iRow++)
// 		for(iCol=0;iCol<nTx*powerProfile.Ndelays();iCol++)
// 		{
// 			tVector sample = LaGenMatDouble::zeros(length,1);
//
// 			// res = mean + L*RandnMatrix(mean.size(),1,0.0,1.0)
// 			Blas_Mat_Vec_Mult(Ls[iCol/nTx],StatUtil::RandnMatrix(length,1,0.0,1.0),sample,1.0,1.0);
//
// 			for(iTime=0;iTime<length;iTime++)
// 				profiledChannelMatrices[iTime](iRow,iCol) = sample(iTime);
// 		}
//
// 	double delay,quotient;
// 	int k;
// 	for(iTime=0;iTime<length;iTime++)
// 	{
// 		for(iRow=0;iRow<nRx;iRow++)
// 			for(iCol=0;iCol<nTx*powerProfile.Ndelays();iCol++)
// 			{
// 				delay = powerProfile.Delays()[iCol/nTx];
// 				quotient = delay / T;
// 				if(quotient == floor(quotient))
// 				{
// #ifdef DEBUG
// 					cout << "se mete en " << int(quotient) << endl;
// #endif
// 					_channelMatrices[iTime](iRow,nTx*int(quotient)+ (iCol % nTx)) += profiledChannelMatrices[iTime](iRow,iCol);
// 					continue;
// 				}
// 				k = int(quotient);
// #ifdef DEBUG
// 				if(k ==1)
// 				{
// 					cout << "iRow = " << iRow << endl;
// 					cout << "delay es " << delay << endl;
// 					cout << "k = " << k << endl;
// 					cout << "profiledChannelMatrices[iTime](iRow,iCol) = " << profiledChannelMatrices[iTime](iRow,iCol) << " iCol = " << iCol << endl;
// 					cout << ((k+1)- quotient)*profiledChannelMatrices[iTime](iRow,iCol) << endl;
// 					cout << (quotient - k)*profiledChannelMatrices[iTime](iRow,iCol) << endl;
// 					cout << "se mete en " << k*nTx + (iCol % nTx) << " y " << (k+1)*nTx + (iCol % nTx) << endl;
// 					cout << "y allí había " << _channelMatrices[iTime](iRow,k*nTx + (iCol % nTx)) << " y " << _channelMatrices[iTime](iRow,(k+1)*nTx + (iCol % nTx)) << endl;
// // 					cout << "Una tecla..."; getchar();
// 				}
// #endif
// 				if(!((k+1) < memory))
// 					throw RuntimeException("BesselChannel::BesselChannel: oooooops...");
// 				_channelMatrices[iTime](iRow,k*nTx + (iCol % nTx)) += ((k+1)- quotient)*profiledChannelMatrices[iTime](iRow,iCol);
// 				_channelMatrices[iTime](iRow,(k+1)*nTx + (iCol % nTx)) += (quotient - k)*profiledChannelMatrices[iTime](iRow,iCol);
//
//
// 			}
// #ifdef DEBUG
// 		cout << "_channelMatrices[iTime] " << endl << _channelMatrices[iTime] << endl;
// 		cout << "profiledChannelMatrices[iTime]" << endl << profiledChannelMatrices[iTime] << endl;
// #endif
// 		_channelMatrices[iTime] = Util::FlipLR(_channelMatrices[iTime]);
// 	}
//
// #ifdef DEBUG
// 	ofstream f("canal",ofstream::trunc);
// 	Util::MatricesVectorToOctaveFileStream(profiledChannelMatrices,"profiledChannelMatrices",f);
// 	Util::ScalarsVectorToOctaveFileStream(powerProfile.Delays(),"retardos",f);
// 	f.close();
// 	cout << "Los retardos" << endl;
// 	Util::Print(powerProfile.Delays());
// #endif
// }

BesselChannel::~BesselChannel()
{
	delete[] _channelMatrices;
}


