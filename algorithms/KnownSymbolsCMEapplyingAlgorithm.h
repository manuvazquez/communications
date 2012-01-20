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
#ifndef KNOWNSYMBOLSCMEAPPLYINGALGORITHM_H
#define KNOWNSYMBOLSCMEAPPLYINGALGORITHM_H

#include <CMEapplyingAlgorithm.h>

/**
	@author Manu <manu@rustneversleeps>
*/
class KnownSymbolsCMEapplyingAlgorithm : public CMEapplyingAlgorithm
{
public:
    KnownSymbolsCMEapplyingAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr,uint N, uint iLastSymbolVectorToBeDetected, vector< ChannelMatrixEstimator * > channelEstimators, MatrixXd preamble, const MatrixXd &symbolVectors);
	virtual bool performsSymbolsDetection() const { return false; }
protected:
    MatrixXd _symbolVectors;


    virtual std::vector<MatrixXd> estimatedChannelMatricesForChannelOrder(uint iChannelOrder,const MatrixXd &observations,const vector<double> &noiseVariances,const MatrixXd& trainingSequence);
    virtual MatrixXd detectedSymbolsForChannelOrder(uint iChannelOrder,const MatrixXd &observations,const vector<double> &noiseVariances,const MatrixXd& trainingSequence);

};

#endif
