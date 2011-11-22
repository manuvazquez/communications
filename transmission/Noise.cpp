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
#include "Noise.h"

Noise::Noise(uint nOutputs,uint length): _nOutputs(nOutputs),_length(length)
{
}

vector<double> Noise::variances() const
{
	vector<double> res(_length);
	for(uint i=0;i<_length;i++)
		res[i] = varianceAt(i);
	return res;
}

MatrixXd Noise::range(uint start,uint end) const
{
	MatrixXd res(_nOutputs,end-start+1);
	
	for(uint i=start;i<=end;i++)
		res.col(i-start) = at(i);
	
	return res;
}
