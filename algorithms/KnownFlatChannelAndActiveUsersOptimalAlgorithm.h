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
#ifndef KNOWNFLATCHANNELANDACTIVEUSERSOPTIMALALGORITHM_H
#define KNOWNFLATCHANNELANDACTIVEUSERSOPTIMALALGORITHM_H

#include <KnownFlatChannelOptimalAlgorithm.h>

/**
It implements the optimal detection algorithm for a flat channel when at each time instant it is perfectly known which users are transmitting (they don't need to be at every moment)

	@author Manu <manu@rustneversleeps>
*/
class KnownFlatChannelAndActiveUsersOptimalAlgorithm : public KnownFlatChannelOptimalAlgorithm
{
public:
    KnownFlatChannelAndActiveUsersOptimalAlgorithm(std::string name, Alphabet alphabet, uint L, uint Nr, uint N, uint iLastSymbolVectorToBeDetected, const MIMOChannel& channel, uint preambleLength, std::vector<std::vector<bool> > usersActivity);

    ~KnownFlatChannelAndActiveUsersOptimalAlgorithm();

protected:
    Alphabet *_noTransmissionAlphabet;
    std::vector<std::vector<bool> > _usersActivity;    
    
    const Alphabet* getAlphabetAt(uint time, int leafHeight) const
    {
        if(_usersActivity[_nInputs-leafHeight][time])
            return &_alphabet;
        else
            return _noTransmissionAlphabet;
    }

};

#endif
