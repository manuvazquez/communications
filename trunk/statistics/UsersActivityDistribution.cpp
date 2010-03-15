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
#include "UsersActivityDistribution.h"

#include <StatUtil.h>

UsersActivityDistribution::UsersActivityDistribution(const double userPersistenceProb, const double newActiveUserProb, const double userPriorProb):_prior(2),_userActiveGivenItWasPdf(2),_userActiveGivenItWasNotPdf(2)
{
    _prior[0] = 1.0 - userPriorProb;
    _prior[1] = userPriorProb;
    
    _userActiveGivenItWasPdf[0] = 1.0 - userPersistenceProb;
    _userActiveGivenItWasPdf[1] = userPersistenceProb;    
    
    _userActiveGivenItWasNotPdf[0] = 1.0 - newActiveUserProb;
    _userActiveGivenItWasNotPdf[1] = newActiveUserProb;    
}

bool UsersActivityDistribution::sampleFromPrior() const
{
    return StatUtil::discrete_rnd(_prior);
}

bool UsersActivityDistribution::sampleGivenItWas(bool previous) const
{
    // the user was active in the last time instant
    if(previous) 
        return bool(StatUtil::discrete_rnd(_userActiveGivenItWasPdf));
    // the user was NOT active in the last time instant
    else
        return bool(StatUtil::discrete_rnd(_userActiveGivenItWasNotPdf));
}

double UsersActivityDistribution::probXgivenY(bool X, bool Y) const
{
    if(X)
        return Y?_userActiveGivenItWasPdf[1]:_userActiveGivenItWasNotPdf[1];
    else
        return Y?_userActiveGivenItWasPdf[0]:_userActiveGivenItWasNotPdf[0];
}

double UsersActivityDistribution::probApriori(bool X) const
{
    if(X)
        return _prior[1];
    else
        return _prior[0];
}

// double UsersActivityDistribution::probApriori(VectorXd &symbolsVector) const
// {
//   double res = 1.0;
//   
//   for(int i=0;i<symbolsVector.size();i++)
// 	res *= probApriori(isUserActive(symbolsVector(i)));
//   
//   return res;
// }
// 
// double UsersActivityDistribution::probXgivenY(VectorXd &X, VectorXd &Y) const
// {
//   if(X.size()!=Y.size())
// 	throw RuntimeException("UsersActivityDistribution::probXgivenY: the sizes of the vector don't match.");
//   
//   double res = 1.0;
//   
//   for(int i=0;i<X.size();i++)
// 	res *= probXgivenY(isUserActive(X(i)),isUserActive(Y(i)));
//   
//   return res;
// }