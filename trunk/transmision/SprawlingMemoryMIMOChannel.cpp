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
#include "SprawlingMemoryMIMOChannel.h"

using namespace std;

SprawlingMemoryMIMOChannel::SprawlingMemoryMIMOChannel(int nTx, int nRx, int length, vector<int> candidateOrders,tMatrix transitionProbabilitiesMatrix,int initialChannelOrderIndex): MIMOChannel(nTx, nRx, length),_candidateOrders(candidateOrders),_transitionProbabilitiesMatrix(transitionProbabilitiesMatrix),_initialChannelOrderIndex(initialChannelOrderIndex)
{
	// it obtains the maximum allowed order for the channel
	int iMaxOrder = 0;
	_maxOrder = _candidateOrders[0];
	for(int i=1;i<_candidateOrders.size();i++)
		if(_candidateOrders[i]>_maxOrder)
		{
			iMaxOrder = i;
			_maxOrder = _candidateOrders[i];
		}

	if(_transitionProbabilitiesMatrix.rows()!=_transitionProbabilitiesMatrix.cols() || _transitionProbabilitiesMatrix.cols()!=_candidateOrders.size())
	{
		throw RuntimeException("SprawlingMemoryMIMOChannel::SprawlingMemoryMIMOChannel: either transition probabilities matrix or channel orders vector dimensions are wrong.");
	}

	// it checks wether the probabilites matrix is coherent
	for(int i=0;i<_transitionProbabilitiesMatrix.rows();i++)
		if(Util::Sum(_transitionProbabilitiesMatrix.row(i))!=1.0)
			throw RuntimeException("SprawlingMemoryMIMOChannel::SprawlingMemoryMIMOChannel: matrix probabilites are not coherent");
}
