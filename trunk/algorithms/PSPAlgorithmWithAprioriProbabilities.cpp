/*
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License version 2 as published by the Free Software Foundation.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public License
   along with this library; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/

#include "PSPAlgorithmWithAprioriProbabilities.h"

void PSPAlgorithmWithAprioriProbabilities::deployState(int iState, const VectorXd& observations, double noiseVariance)
{
PSPAlgorithm::deployState(iState, observations, noiseVariance);
}

void PSPAlgorithmWithAprioriProbabilities::run(MatrixXd observations, vector< double > noiseVariances)
{
//     if(observations.cols()<(_startDetectionTime+1+_d))
//         throw RuntimeException("PSPAlgorithm::Run: not enough observations.");
// 
//     // the last N*(m-1) symbols of the preamble are copied into a c++ vector...
//     int preambleLength = _preamble.rows()*_preamble.cols();
//     vector<tSymbol> initialStateVector(_nInputs*(_channelOrder-1));
// 
//     // (it must be taken into account that the number of columns of the preamble might be greater than m-1)
//     int iFirstPreambleSymbolNeeded = (_preamble.cols()-(_channelOrder-1))*_nInputs;
//     for(int i=iFirstPreambleSymbolNeeded;i<preambleLength;i++)
//         initialStateVector[i-iFirstPreambleSymbolNeeded] = _preamble(i % _nInputs,i / _nInputs);
// 
//     // ...in order to use the method "SymbolsVectorToInt" from "Alphabet" to obtain the initial state
//     int initialState = _alphabet.symbolsArray2int(initialStateVector);
// 
//     // the initial state is initalized
//     _exitStage[initialState][0] = PSPPath(_iLastSymbolVectorToBeDetected+_d,0.0,_preamble,vector<vector<MatrixXd> > (1,vector<MatrixXd>(0)),vector<ChannelMatrixEstimator *>(1,_channelEstimator));
// 
//     process(observations,noiseVariances);
}

void PSPAlgorithmWithAprioriProbabilities::run(MatrixXd observations, vector< double > noiseVariances, MatrixXd trainingSequence)
{
  throw RuntimeException("PSPAlgorithmWithAprioriProbabilities::run: this is not implemented.");
}

