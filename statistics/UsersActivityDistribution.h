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
#ifndef USERSACTIVITYDISTRIBUTION_H
#define USERSACTIVITYDISTRIBUTION_H

/**
It implements...

	@author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <vector>

class UsersActivityDistribution{
protected:
    std::vector<double> _prior;
    std::vector<double> _userActiveGivenItWasPdf;
    std::vector<double> _userActiveGivenItWasNotPdf;
public:
    UsersActivityDistribution(const double userPersistenceProb, const double newActiveUserProb, const double userPriorProb);
    bool sampleFromPrior() const;
    bool sampleGivenItWas(bool previous) const;
    double probXgivenY(bool X, bool Y) const;
    double probApriori(bool X) const;
	void setApriori(const double userPriorProb) { _prior[0] = 1.0 - userPriorProb; _prior[1] = userPriorProb;}
};

#endif
