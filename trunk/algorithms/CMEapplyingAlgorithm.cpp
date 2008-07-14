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
#include "CMEapplyingAlgorithm.h"

// #define EXPORT_REAL_DATA

#ifdef EXPORT_REAL_DATA
    extern MIMOChannel *realChannel;
    extern tMatrix *realSymbols;
    extern Noise *realNoise;
#endif

CMEapplyingAlgorithm::CMEapplyingAlgorithm(string name, Alphabet alphabet, int L, int N, int K, vector< ChannelMatrixEstimator * > channelEstimators, tMatrix preamble): UnknownChannelOrderAlgorithm(name, alphabet, L, N, K, channelEstimators, preamble, preamble.cols())
{
}

void CMEapplyingAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    throw RuntimeException("CMEapplyingAlgorithm::Run (without training sequence) not implemented.");
}

void CMEapplyingAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
    int m,iTxAntenna,iRxAntenna,iDelay;
    tRange rAll;
    tVector CMEs(_candidateOrders.size());

#ifdef EXPORT_REAL_DATA
    tMatrix channelMatrix = (*realChannel)[_preamble.cols()];
#endif

    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
    {
        vector<tMatrix> estimatedChannelMatrices = estimatedChannelMatricesForChannelOrder(iChannelOrder,observations,noiseVariances,trainingSequence);
        tMatrix detectedSymbolVectors = detectedSymbolsForChannelOrder(iChannelOrder,observations,noiseVariances,trainingSequence);

        tMatrix preambleDetectedSymbolVectors = Util::Append(_preamble,detectedSymbolVectors);

        int nSymbolVectors = detectedSymbolVectors.cols();
        double variance = noiseVariances[detectedSymbolVectors.cols()-1];

        m = _candidateOrders[iChannelOrder];

        tMatrix estimatedChannelMatrix = estimatedChannelMatrices[estimatedChannelMatrices.size()-1];

        vector<tVector> hs(_L,LaGenMatDouble::zeros(_N*m,1));

        tMatrix C(nSymbolVectors,_N*m);
        for(iTxAntenna=0;iTxAntenna<_N;iTxAntenna++)
            for(iDelay=0;iDelay<m;iDelay++)
            {
                // symbols are transformed
                for(int CmatrixRow=0;CmatrixRow<nSymbolVectors;CmatrixRow++)
                    C(CmatrixRow,iTxAntenna*m+iDelay) = preambleDetectedSymbolVectors(iTxAntenna,_preamble.cols()-iDelay+CmatrixRow);

                // channel is transformed
                for(iRxAntenna=0;iRxAntenna<_L;iRxAntenna++)
                    hs[iRxAntenna](iTxAntenna*m+iDelay) = estimatedChannelMatrix(iRxAntenna,iTxAntenna+(m-1-iDelay)*_N);
            }

        // CME
        double CME = 0.0;
        tRange rAllObservationsCols(_preamble.cols(),preambleDetectedSymbolVectors.cols()-1);

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

tMatrix CMEapplyingAlgorithm::GetDetectedSymbolVectors()
{
    return tMatrix(0,0);
}

vector<tMatrix> CMEapplyingAlgorithm::GetEstimatedChannelMatrices()
{
    return vector<tMatrix>(0);
}

