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
    extern MatrixXd *realSymbols;
    extern Noise *realNoise;
#endif

CMEapplyingAlgorithm::CMEapplyingAlgorithm(string name, Alphabet alphabet, int L, int Nr,int N, int iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, MatrixXd preamble): UnknownChannelOrderAlgorithm(name, alphabet, L, Nr,N, iLastSymbolVectorToBeDetected, channelEstimators, preamble, preamble.cols())
{
}

void CMEapplyingAlgorithm::run(MatrixXd observations,vector<double> noiseVariances)
{
    throw RuntimeException("CMEapplyingAlgorithm::run (without training sequence) not implemented.");
}

void CMEapplyingAlgorithm::run(MatrixXd observations,vector<double> noiseVariances, MatrixXd trainingSequence)
{
    int m,iTxAntenna,iRxAntenna,iDelay;
    tVector CMEs(_candidateOrders.size());

#ifdef EXPORT_REAL_DATA
    MatrixXd channelMatrix = realChannel->at(_preamble.cols());
#endif

    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
    {
        vector<MatrixXd> estimatedChannelMatrices = estimatedChannelMatricesForChannelOrder(iChannelOrder,observations,noiseVariances,trainingSequence);
        MatrixXd detectedSymbolVectors = detectedSymbolsForChannelOrder(iChannelOrder,observations,noiseVariances,trainingSequence);

        MatrixXd preambleDetectedSymbolVectors(_preamble.rows(),_preamble.cols()+detectedSymbolVectors.cols());
        
//         preambleDetectedSymbolVectors << Util::lapack2eigen(_preamble),detectedSymbolVectors;
        preambleDetectedSymbolVectors << _preamble,detectedSymbolVectors;

        int nSymbolVectors = detectedSymbolVectors.cols();
        double variance = noiseVariances[detectedSymbolVectors.cols()-1];

        m = _candidateOrders[iChannelOrder];

        MatrixXd estimatedChannelMatrix = estimatedChannelMatrices[estimatedChannelMatrices.size()-1];

        vector<VectorXd> hs(_nOutputs,VectorXd::Zero(_nInputs*m,1));

        MatrixXd C(nSymbolVectors,_nInputs*m);
        for(iTxAntenna=0;iTxAntenna<_nInputs;iTxAntenna++)
            for(iDelay=0;iDelay<m;iDelay++)
            {
                // symbols are transformed
                for(int CmatrixRow=0;CmatrixRow<nSymbolVectors;CmatrixRow++)
                    C(CmatrixRow,iTxAntenna*m+iDelay) = preambleDetectedSymbolVectors(iTxAntenna,_preamble.cols()-iDelay+CmatrixRow);

                // channel is transformed
                for(iRxAntenna=0;iRxAntenna<_nOutputs;iRxAntenna++)
                    hs[iRxAntenna](iTxAntenna*m+iDelay) = estimatedChannelMatrix(iRxAntenna,iTxAntenna+(m-1-iDelay)*_nInputs);
            }

        // CME
        double CME = 0.0;

        for(iRxAntenna=0;iRxAntenna<_nOutputs;iRxAntenna++)
        {
            VectorXd error = observations.block(iRxAntenna,_preamble.cols(),1,preambleDetectedSymbolVectors.cols()-_preamble.cols()) - C*hs[iRxAntenna];

            // CME += error'*error
            CME += error.dot(error);
        }
        CME /= variance;

        CME += _nOutputs*log(fabs((C.transpose()*C).determinant()));
        CME /= 2.0;

        CMEs(iChannelOrder) = CME;
    } // for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)

    for(uint iChannelOrder=0;iChannelOrder<_candidateOrders.size();iChannelOrder++)
        _channelOrderAPPs.row(iChannelOrder).setConstant(CMEs(iChannelOrder));
}

MatrixXd CMEapplyingAlgorithm::getDetectedSymbolVectors_eigen()
{
    return MatrixXd(0,0);
}

vector<MatrixXd> CMEapplyingAlgorithm::getEstimatedChannelMatrices_eigen()
{
    return vector<MatrixXd>(0);
}

