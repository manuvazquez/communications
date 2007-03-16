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
#ifndef PSPPATH_H
#define PSPPATH_H

#include <ViterbiPath.h>

/**
	@author Manu <manu@rustneversleeps>
*/

// #define DEBUG13
#include <vector>
#include <ChannelMatrixEstimator.h>

class PSPPath : public ViterbiPath
{
protected:
	std::vector<ChannelMatrixEstimator *> _channelMatrixEstimators;
	tMatrix **_estimatedChannelMatrices;
public:
    PSPPath();

    PSPPath(int nTimeInstants,double cost, tMatrix initialSequence, std::vector<std::vector<tMatrix> > initialChannelMatrices, std::vector<ChannelMatrixEstimator *> channelMatrixEstimators);

    PSPPath(const PSPPath &path);

    ~PSPPath();

	ChannelMatrixEstimator * GetChannelMatrixEstimator() const { return _channelMatrixEstimators[0];}
	tMatrix GetChannelMatrix(int n)
	{
		#ifdef DEBUG13
			cout << "nº de vectores de símbolos detectados " << _detectedSequence->cols() << endl;
			cout << "DEvolviendo" << endl << _estimatedChannelMatrices[0][n] << endl;
		#endif
		return _estimatedChannelMatrices[0][n];
	}
    void Clean();
    void Print() const;
    /**
     * Updates the current path object from another path object, plus a new symbol vector, a new cost, a new vector of ChannelMatrixEstimators, and a vector of channel matrices
     * @param path
     * @param newSymbolVector
     * @param newCost
     * @param newChannelMatrixEstimators the estimators are directly stored (they are not cloned)
     * @param newChannelMatrices
     */
    void Update(const PSPPath& path, tVector newSymbolVector, double newCost, std::vector<ChannelMatrixEstimator *> newChannelMatrixEstimators/*,const std::vector<tMatrix> &newChannelMatrices*/);
	void operator=(const PSPPath &path);

};

#endif
