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
#ifndef KNOWNCHANNELALGORITHM_H
#define KNOWNCHANNELALGORITHM_H

#include <Algorithm.h>
#include <MIMOChannel.h>

/**
	@author Manu <manu@rustneversleeps>
*/

class KnownChannelAlgorithm : public Algorithm
{
protected:
	const MIMOChannel &_channel;
public:
    KnownChannelAlgorithm(string name, Alphabet alphabet,int L,int Nr,int N, int iLastSymbolVectorToBeDetected, const MIMOChannel &channel);
    
    using Algorithm::run;

    void run(tMatrix observations,vector<double> noiseVariances, tMatrix trainingSequence);
    vector<tMatrix> getEstimatedChannelMatrices();
};

#endif
