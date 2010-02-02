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
#ifndef CHANNELDEPENDENTNOISE_H
#define CHANNELDEPENDENTNOISE_H

#include <Noise.h>
#include <math.h>
#include <StillMemoryMIMOChannel.h>

/**
    @author Manu <manu@rustneversleeps>
*/
class ChannelDependentNoise : public Noise
{
protected:
    MatrixXd _matrix;
    MIMOChannel *_channel;
    std::vector<double> _stdDevs;
public:
    ChannelDependentNoise(MIMOChannel *channel);

    virtual void setSNR(int SNR,double alphabetVariance);
    virtual void print() const { cout << _matrix;}
    double stdDevAt(int n) const;
    VectorXd at(uint n) const;
    virtual MatrixXd range_eigen(int start,int end) const;
};

#endif