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
#include "CMEBasedAlgorithm.h"

// #define EXPORT_REAL_DATA

#ifdef EXPORT_REAL_DATA
	extern MIMOChannel *realChannel;
	extern tMatrix *realSymbols;
	extern Noise *realNoise;
#endif

CMEBasedAlgorithm::CMEBasedAlgorithm(string name, Alphabet alphabet, int L, int N, int frameLength, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, const tMatrix &symbolVectors): UnknownChannelOrderAlgorithm(name, alphabet, L, N, frameLength, channelEstimators, preamble, iFirstObservation),_symbolVectors(symbolVectors)
{
}

void CMEBasedAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
	int m,iTxAntenna,iRxAntenna,iDelay;
	int nSymbolVectors = _symbolVectors.cols() - _preamble.cols();
	tRange rAll;
	tVector CMEs(_candidateOrders.size());

#ifdef EXPORT_REAL_DATA
	tMatrix channelMatrix = (*realChannel)[_preamble.cols()];
#endif
	double variance = noiseVariances[_symbolVectors.cols()-1];

	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
	{
		m = _candidateOrders[iChannelOrder];

		// channel estimation
		tRange rSymbolVectors(_preamble.cols()-m+1,_preamble.cols());
		for(int iSymbolVector=_preamble.cols();iSymbolVector<_K;iSymbolVector++)
		{
			_channelEstimators[iChannelOrder]->NextMatrix(observations.col(iSymbolVector),_symbolVectors(rAll,rSymbolVectors),noiseVariances[iSymbolVector]);
			rSymbolVectors = rSymbolVectors + 1;
		}

		tMatrix estimatedChannelMatrix = _channelEstimators[iChannelOrder]->LastEstimatedChannelMatrix();

		vector<tVector> hs(_L,LaGenMatDouble::zeros(_N*m,1));

		tMatrix C(nSymbolVectors,_N*m);
		for(iTxAntenna=0;iTxAntenna<_N;iTxAntenna++)
			for(iDelay=0;iDelay<m;iDelay++)
			{
				// symbols are transformed
				for(int CmatrixRow=0;CmatrixRow<nSymbolVectors;CmatrixRow++)
					C(CmatrixRow,iTxAntenna*m+iDelay) = _symbolVectors(iTxAntenna,_preamble.cols()-iDelay+CmatrixRow);

				// channel is transformed
				for(iRxAntenna=0;iRxAntenna<_L;iRxAntenna++)
					hs[iRxAntenna](iTxAntenna*m+iDelay) = estimatedChannelMatrix(iRxAntenna,iTxAntenna+(m-1-iDelay)*_N);
			}

		// CME
		double CME = 0.0;
		tRange rAllObservationsCols(_preamble.cols(),_symbolVectors.cols()-1);
		for(iRxAntenna=0;iRxAntenna<_L;iRxAntenna++)
		{
			// error = R
			tVector error = observations(tRange(iRxAntenna),rAllObservationsCols);

			// error = error - C * hs[iRxAntenna]
			Blas_Mat_Vec_Mult(C,hs[iRxAntenna],error,-1.0,1.0);

			// CME += error'*error
			CME += Blas_Dot_Prod(error,error);
		}
		CME /= variance;

		tMatrix CTransC(_N*m,_N*m);

		//  CTransC = C'*C
		Blas_Mat_Trans_Mat_Mult(C,C,CTransC);

		// LU decomposition is applied: in CTransC wil now be U
		tLongIntVector piv(_N*m);
		LUFactorizeIP(CTransC,piv);

		double detCTransC = 1.0;
		for(int iDiag=0;iDiag<CTransC.cols();iDiag++)
			detCTransC *= CTransC(iDiag,iDiag);

		CME += _L*log(fabs(detCTransC));
		CME /= 2.0;

		CMEs(iChannelOrder) = CME;
	} // for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)

	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
		_channelOrderAPPs.row(iChannelOrder) = CMEs(iChannelOrder);
}

void CMEBasedAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
	Run(observations,noiseVariances);
}

tMatrix CMEBasedAlgorithm::getDetectedSymbolVectors()
{
	return tMatrix(0,0);
}

vector<tMatrix> CMEBasedAlgorithm::GetEstimatedChannelMatrices()
{
	return vector<tMatrix>(0);
}

