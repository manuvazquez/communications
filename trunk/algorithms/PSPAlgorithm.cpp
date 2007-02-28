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
#include "PSPAlgorithm.h"

#define DEBUG

PSPAlgorithm::PSPAlgorithm(string name, Alphabet alphabet, int L, int N, int K, int m, ChannelMatrixEstimator* channelEstimator, tMatrix preamble, int smoothingLag): KnownChannelOrderAlgorithm(name, alphabet, L, N, K, m, channelEstimator, preamble),_d(smoothingLag),_startDetectionTime(preamble.cols()),_trellis(alphabet,N,m),_detectedSymbolVectors(NULL)
{
    if(preamble.cols() < (m-1))
        throw RuntimeException("PSPAlgorithm::PSPAlgorithm: preamble dimensions are wrong.");

    _exitStage = new PSPPath[_trellis.Nstates()];
    _arrivalStage = new PSPPath[_trellis.Nstates()];
}


PSPAlgorithm::~PSPAlgorithm()
{
    delete[] _exitStage;
    delete[] _arrivalStage;
    delete _detectedSymbolVectors;
}

void PSPAlgorithm::Run(tMatrix observations,vector<double> noiseVariances)
{
    if(observations.cols()<(_startDetectionTime+1+_d))
        throw RuntimeException("PSPAlgorithm::Run: Not enough observations.");
}

void PSPAlgorithm::Run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence)
{
    if(observations.rows()!=_L || trainingSequence.rows()!=_N)
        throw RuntimeException("PSPAlgorithm::Run: Observations matrix or training sequence dimensions are wrong.");

    // to process the training sequence, we need both the preamble and the symbol vectors related to it
    tMatrix preambleTrainingSequence = Util::Append(_preamble,trainingSequence);


    tRange rSymbolVectorsTrainingSequece(0,preambleTrainingSequence.cols()-1);

    vector<tMatrix> trainingSequenceChannelMatrices = ProcessTrainingSequence(observations,noiseVariances,trainingSequence);

	#ifdef DEBUG
    	cout << "cols es " << preambleTrainingSequence.cols() << endl;
   	#endif

   	PSPPath prueba(_K+_d,0.0,preambleTrainingSequence,vector<vector<tMatrix> > (1,trainingSequenceChannelMatrices),vector<ChannelMatrixEstimator *>(1,_channelEstimator));

   	prueba.Print();

    // the symbols contained in the preamble are copied into a c++ vector...
    int preambleTrainingSequenceLength = preambleTrainingSequence.rows()*preambleTrainingSequence.cols();
    vector<tSymbol> initialStateVector(_N*(_m-1));

    // (it must be taken into account that the number of columns of the preamble might be greater than m-1)
    int iFirstPreambleSymbolNeeded = (preambleTrainingSequence.cols()-(_m-1))*_N;
    for(int i=iFirstPreambleSymbolNeeded;i<preambleTrainingSequenceLength;i++)
        initialStateVector[i-iFirstPreambleSymbolNeeded] = preambleTrainingSequence(i % _N,i / _N);

    // ...in order to use the method "SymbolsVectorToInt" from "Alphabet"
    int initialState = _alphabet.SymbolsArrayToInt(initialStateVector);

	#ifdef DEBUG
    	cout << "El estado inicial es " << initialState << endl;
   	#endif
}

void PSPAlgorithm::Process(const tMatrix &observations,vector<double> noiseVariances)
{
    // memory for the symbol vectors being detected is reserved
    _detectedSymbolVectors = new tMatrix(_N,_K+_d);
}

tMatrix PSPAlgorithm::GetDetectedSymbolVectors()
{
	return tMatrix(0,0);
}

vector<tMatrix> PSPAlgorithm::GetEstimatedChannelMatrices()
{
	return vector<tMatrix>(1,tMatrix(0,0));
}
