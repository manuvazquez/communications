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
#include "StatUtil.h"

vector<int> StatUtil::Discrete_rnd(int nSamples, tVector probabilities,Random &randomGenerator)
{
    int i,j;
	double uniform;

    tVector normalizedProbabilities = Util::Normalize(probabilities);
    int nProbabilities = probabilities.size();
    
    double *distributionFunction = new double[nProbabilities];
    double acum = 0.0;
    for(i=0;i<nProbabilities;i++)
    {
        acum += normalizedProbabilities(i);
        distributionFunction[i] = acum;
    }

    vector<int> res(nSamples);

    for(i=0;i<nSamples;i++)
    {
		uniform = randomGenerator.rand();
		j=0;
		while(uniform>distributionFunction[j])
			j++;
		res[i] = j;
    }
	
	// memory release
	delete[] distributionFunction;
	return res;
}

tMatrix StatUtil::RandnMatrix(int rows,int cols,double mean,double variance,Random &randomGenerator)
{
	tMatrix res(rows,cols);
	double stdDv = sqrt(variance);

	int j;
	for(int i=0;i<rows;i++)
		for(j=0;j<cols;j++)
			res(i,j) = 	randomGenerator.randn()*stdDv + mean;

	return res;
}
