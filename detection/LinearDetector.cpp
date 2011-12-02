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
#include "LinearDetector.h"

LinearDetector::LinearDetector(uint rows,uint cols,double alphabetVariance):_channelMatrixRows(rows),_channelMatrixCols(cols),_alphabetVariance(alphabetVariance)
{
}

void LinearDetector::stateStepsFromObservationsSequence(const MatrixXd &observations,uint smoothingLag,uint iFrom,uint iTo)
{
    for(uint i=iFrom;i<iTo;i++)
        stateStep(Util::toVector(observations.block(0,i,observations.rows(),smoothingLag+1),columnwise));
}