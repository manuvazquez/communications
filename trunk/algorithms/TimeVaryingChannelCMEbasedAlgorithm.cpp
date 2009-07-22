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
#include "TimeVaryingChannelCMEbasedAlgorithm.h"

// #define EXPORT_REAL_DATA

#ifdef EXPORT_REAL_DATA
	extern MIMOChannel *realChannel;
	extern tMatrix *realSymbols;
	extern Noise *realNoise;
#endif

TimeVaryingChannelCMEbasedAlgorithm::TimeVaryingChannelCMEbasedAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, const tMatrix &symbolVectors): UnknownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, preamble, iFirstObservation),_symbolVectors(symbolVectors)
{
}

void TimeVaryingChannelCMEbasedAlgorithm::run(tMatrix observations,vector<double> noiseVariances)
{
	int m,iTxAntenna,iDelay;
	int nSymbolVectors = _symbolVectors.cols() - _preamble.cols();
	tRange rAll;
	tVector CMEs(_candidateOrders.size());
	tVector observationError;
	double accumulatedSquaredObservationsError;

#ifdef EXPORT_REAL_DATA
	tMatrix channelMatrix = (*realChannel)[_preamble.cols()];
#endif

	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
	{
		m = _candidateOrders[iChannelOrder];

		accumulatedSquaredObservationsError = 0.0;

		// channel estimation
		tRange rSymbolVectors(_preamble.cols()-m+1,_preamble.cols());
		int skipNumber = 0;
		for(int iSymbolVector=_preamble.cols();iSymbolVector<observations.cols();iSymbolVector++)
		{
			_channelEstimators[iChannelOrder]->nextMatrix(observations.col(iSymbolVector),_symbolVectors(rAll,rSymbolVectors),noiseVariances[iSymbolVector]);

			observationError = observations.col(iSymbolVector);
			// observationError = observationError - _channelEstimators[iChannelOrder]->lastEstimatedChannelMatrix() * _symbolVectors(rAll,rSymbolVectors)
			Blas_Mat_Vec_Mult(_channelEstimators[iChannelOrder]->lastEstimatedChannelMatrix(),Util::toVector(_symbolVectors(rAll,rSymbolVectors),columnwise),observationError,-1.0,1.0);

			accumulatedSquaredObservationsError += double(skipNumber>50)*Blas_Dot_Prod(observationError,observationError)/noiseVariances[iSymbolVector];

			skipNumber++;

			rSymbolVectors = rSymbolVectors + 1;
		}

		tMatrix estimatedChannelMatrix = _channelEstimators[iChannelOrder]->lastEstimatedChannelMatrix();

		vector<tVector> hs(_nOutputs,LaGenMatDouble::zeros(_nInputs*m,1));

		tMatrix C(nSymbolVectors,_nInputs*m);
		for(iTxAntenna=0;iTxAntenna<_nInputs;iTxAntenna++)
			for(iDelay=0;iDelay<m;iDelay++)
				// symbols are transformed
				for(int CmatrixRow=0;CmatrixRow<nSymbolVectors;CmatrixRow++)
					C(CmatrixRow,iTxAntenna*m+iDelay) = _symbolVectors(iTxAntenna,_preamble.cols()-iDelay+CmatrixRow);

		// CME
		double CME = 0.0;

		CME = accumulatedSquaredObservationsError;

		tMatrix CTransC(_nInputs*m,_nInputs*m);

		//  CTransC = C'*C
		Blas_Mat_Trans_Mat_Mult(C,C,CTransC);

		// LU decomposition is applied: in CTransC wil now be U
		tLongIntVector piv(_nInputs*m);
		LUFactorizeIP(CTransC,piv);

		double detCTransC = 1.0;
		for(int iDiag=0;iDiag<CTransC.cols();iDiag++)
			detCTransC *= CTransC(iDiag,iDiag);

		CME += _nOutputs*log(fabs(detCTransC));
		CME /= 2.0;

		CMEs(iChannelOrder) = CME;
	} // for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)

	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
		_channelOrderAPPs.row(iChannelOrder) = CMEs(iChannelOrder);
}

void TimeVaryingChannelCMEbasedAlgorithm::run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
	run(observations,noiseVariances);
}

tMatrix TimeVaryingChannelCMEbasedAlgorithm::getDetectedSymbolVectors()
{
	return tMatrix(0,0);
}

vector<tMatrix> TimeVaryingChannelCMEbasedAlgorithm::getEstimatedChannelMatrices()
{
	return vector<tMatrix>(0);
}

