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
#ifndef ARMULTIUSERCDMACHANNEL_H
#define ARMULTIUSERCDMACHANNEL_H

#include <MultiuserCDMAchannel.h>
#include <ARprocess.h>

/**
It models a multiuser CDMA channel in which each user has a single coefficent (tantamount to transmitting power). The different coefficients evolve according to a order 2 AutoRegressive process.

    @author Manu <manu@rustneversleeps>
*/
class ARMultiuserCDMAchannel : public MultiuserCDMAchannel
{
public:
    ARMultiuserCDMAchannel(int length, const tMatrix& spreadingCodes,const ARprocess &arProcess);

protected:
    ARprocess _ARprocess;
    vector<MatrixXd> _userCoeffs;
    
    virtual MatrixXd at(int n) const { return _userCoeffs[n];}
};

#endif
