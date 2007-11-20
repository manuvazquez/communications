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

// #define DEBUG5
// #define DEBUG7
// #define EXPORT_REAL_DATA

#ifdef EXPORT_REAL_DATA
	extern MIMOChannel *realChannel;
	extern tMatrix *realSymbols;
	extern Noise *realNoise;
#endif

TimeVaryingChannelCMEbasedAlgorithm::TimeVaryingChannelCMEbasedAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble, int iFirstObservation, const tMatrix &symbolVectors): UnknownChannelOrderAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, iFirstObservation),_symbolVectors(symbolVectors)
{
}

void TimeVaryingChannelCMEbasedAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
	int m,iTxAntenna,iRxAntenna,iDelay;
	int nSymbolVectors = _symbolVectors.cols() - _preamble.cols();
	tRange rAll;
	tVector CMEs(_candidateOrders.size());
	tVector observationError;
	double accumulatedSquaredObservationsError;

#ifdef DEBUG2
	cout << "nº = " << nSymbolVectors << endl;
	cout << "el primer canal es (sigma2 = " << noiseVariances[_preamble.cols()] << ")" << endl << (*realChannel)[_preamble.cols()];
	cout << "el ultimo canal es (sigma2 = " << noiseVariances[_symbolVectors.cols()-1] << ")" << endl << (*realChannel)[_symbolVectors.cols()-1];
	cout << "Una tecla..."; getchar();
#endif

#ifdef EXPORT_REAL_DATA
	tMatrix channelMatrix = (*realChannel)[_preamble.cols()];
#endif
	double variance = noiseVariances[_symbolVectors.cols()-1];

#ifdef DEBUG3
	for(int hola=_preamble.cols();hola<_K;hola++)
		cout << "varianza en " << hola << " = " << noiseVariances[hola] << endl;
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
			_channelEstimators[iChannelOrder]->NextMatrix(observations.col(iSymbolVector),_symbolVectors(rAll,rSymbolVectors),noiseVariances[iSymbolVector]);

			observationError = observations.col(iSymbolVector);
			// observationError = observationError - _channelEstimators[iChannelOrder]->LastEstimatedChannelMatrix() * _symbolVectors(rAll,rSymbolVectors)
			Blas_Mat_Vec_Mult(_channelEstimators[iChannelOrder]->LastEstimatedChannelMatrix(),Util::ToVector(_symbolVectors(rAll,rSymbolVectors),columnwise),observationError,-1.0,1.0);

			accumulatedSquaredObservationsError += double(skipNumber>50)*Blas_Dot_Prod(observationError,observationError)/noiseVariances[iSymbolVector];

			skipNumber++;

			rSymbolVectors = rSymbolVectors + 1;
		}

		tMatrix estimatedChannelMatrix = _channelEstimators[iChannelOrder]->LastEstimatedChannelMatrix();

// 		rSymbolVectors.set(_preamble.cols()-m+1,_preamble.cols(),1);
// 		for(int iSymbolVector=_preamble.cols();iSymbolVector<observations.cols();iSymbolVector++)
// 		{
// 			observationError = observations.col(iSymbolVector);
// 			// observationError = observationError - _channelEstimators[iChannelOrder]->LastEstimatedChannelMatrix() * _symbolVectors(rAll,rSymbolVectors)
// 			Blas_Mat_Vec_Mult(estimatedChannelMatrix,Util::ToVector(_symbolVectors(rAll,rSymbolVectors),columnwise),observationError,-1.0,1.0);
//
// #ifdef DEBUG5
// 			tVector observationWithouNoise = observationError;
// 			observationWithouNoise = 0.0;
// 			Blas_Mat_Vec_Mult(estimatedChannelMatrix,Util::ToVector(_symbolVectors(rAll,rSymbolVectors),columnwise),observationWithouNoise);
// 			cout << "observationWithouNoise = " << endl << observationWithouNoise;
// #endif
//
// 			accumulatedSquaredObservationsError += Blas_Dot_Prod(observationError,observationError);
//
// 			rSymbolVectors = rSymbolVectors + 1;
// 		}

#ifdef DEBUG
		cout << "------------- m = " << m << " --------------------" << endl;
		cout << "estimatedChannelMatrix es" << endl << estimatedChannelMatrix;
#endif

#ifdef DEBUG2
		cout << "estimatedChannelMatrix es" << endl << estimatedChannelMatrix;
		cout << "el canal de verdad es" << endl << channelMatrix;
		cout << "Una tecla..."; getchar();
#endif

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


// 		tRange rAllObservationsCols(_preamble.cols(),_symbolVectors.cols()-1);
// 		for(iRxAntenna=0;iRxAntenna<_L;iRxAntenna++)
// 		{
// 			// error = R
// 			tVector error = observations(tRange(iRxAntenna),rAllObservationsCols);
//
// 			// error = error - C * hs[iRxAntenna]
// 			Blas_Mat_Vec_Mult(C,hs[iRxAntenna],error,-1.0,1.0);
// #ifdef DEBUG2
// 			cout << "error es" << endl << error;
// #endif
//
// #ifdef DEBUG5
// 			tVector Ch(_symbolVectors.cols()-_preamble.cols());
// 			Blas_Mat_Vec_Mult(C,hs[iRxAntenna],Ch);
// 			cout << "Ch es" << endl << Ch;
// 			cout << "Una tecla..."; getchar();
// #endif
//
// #ifdef DEBUG4
// 			cout << "el tamaño de R es " << error.size() << endl;
// 			cout << "nSymbolVectors = " << nSymbolVectors << endl;
// #endif
// 			// CME += error'*error
// 			CME += Blas_Dot_Prod(error,error);
// 		}

#ifdef DEBUG6
		cout << "CME = " << CME << " y accumulatedSquaredObservationsError = " << accumulatedSquaredObservationsError << endl;
		cout << "CME==accumulatedSquaredObservationsError: " << (CME == accumulatedSquaredObservationsError) << endl;
#endif



		CME = accumulatedSquaredObservationsError;

// 		CME /= variance;

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

#ifdef DEBUG7
		cout << "El CME calculado es " << CME << endl;
#endif

#ifdef DEBUG2
		cout << "la varianza es " << variance << endl;
		cout << "CME (m = " << m << ")" << " = " << endl << CME << endl;
		cout << "Una tecla..."; getchar();
#endif
	} // for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)

#ifdef DEBUG_PRINT_CMES
	cout << "Los CMEs son" << endl << CMEs;
	cout << "y normalizados.." << endl << normalizedInvCMEs;
#endif

	for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
		_channelOrderAPPs.row(iChannelOrder) = CMEs(iChannelOrder);
}

void TimeVaryingChannelCMEbasedAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
	Run(observations,noiseVariances);
}

tMatrix TimeVaryingChannelCMEbasedAlgorithm::GetDetectedSymbolVectors()
{
	return tMatrix(0,0);
}

vector<tMatrix> TimeVaryingChannelCMEbasedAlgorithm::GetEstimatedChannelMatrices()
{
	return vector<tMatrix>(0);
}

