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

// the seed used to create the random objects is generated from the system time
#define RANDOM_SEED

using namespace std;

int StatUtil::Discrete_rnd(const tVector &probabilities)
{
    int i;
	double uniform;

	#ifdef RANDOM_SEED
		static Random randomGenerator;
	#else
		static Random randomGenerator(11);
	#endif

//     tVector normalizedProbabilities = Util::Normalize(probabilities);
    int nProbabilities = probabilities.size();

    double *distributionFunction = new double[nProbabilities];
    distributionFunction[0] = probabilities(0);
    for(i=1;i<nProbabilities;i++)
           distributionFunction[i] = distributionFunction[i-1]+probabilities(i);

	uniform = randomGenerator.rand();
	int res = 0;
	while(uniform>distributionFunction[res])
		res++;

	// memory release
	delete[] distributionFunction;

	return res;
}

vector<int> StatUtil::Discrete_rnd(int nSamples,const tVector &probabilities)
{
    int i,j;
	double uniform;

	#ifdef RANDOM_SEED
		static Random randomGenerator;
	#else
		static Random randomGenerator(1);
	#endif

    tVector normalizedProbabilities = Util::Normalize(probabilities);
    int nProbabilities = probabilities.size();

    double *distributionFunction = new double[nProbabilities];
    distributionFunction[0] = normalizedProbabilities(0);
    for(i=1;i<nProbabilities;i++)
           distributionFunction[i] = distributionFunction[i-1]+normalizedProbabilities(i);

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

tMatrix StatUtil::RandnMatrix(int rows,int cols,double mean,double variance)
{
	tMatrix res(rows,cols);
	double stdDv = sqrt(variance);

	#ifdef RANDOM_SEED
		static Random randomGenerator;
	#else
		static Random randomGenerator(2);
	#endif

	int j;
	for(int i=0;i<rows;i++)
		for(j=0;j<cols;j++)
			res(i,j) = 	randomGenerator.randn()*stdDv + mean;

	return res;
}

double StatUtil::NormalPdf(double x,double mean,double variance)
{
	double distance = fabs(x - mean);

	return 1.0/sqrt(2.0*M_PI*variance)*exp(-(distance*distance)/(2.0*variance));
}

double StatUtil::NormalPdf(const tVector &x,const tVector &mean,const tMatrix &covariance)
{

	int N = x.size();
	// the received covariance matrix can't be modified
	tMatrix invCovariance = covariance;

	tLongIntVector piv(N);
	LUFactorizeIP(invCovariance,piv);
	double detCovariance = 1.0;
	for(int i=0;i<N;i++)
		detCovariance *= invCovariance(i,i);

	// invCovariance = inv(covariance)
	LaLUInverseIP(invCovariance,piv);

	tVector xMinusMean(N);
	// xMinusMean = x - mean
	Util::Add(x,mean,xMinusMean,1.0,-1.0);

	tVector invCovarianceXminusMean(N);
	// invCovarianceXminusMean = -0.5 * invCovariance * xMinusMean
	Blas_Mat_Vec_Mult(invCovariance,xMinusMean,invCovarianceXminusMean,-0.5);

	return 1.0/(sqrt(fabs(detCovariance))*pow(2.0*M_PI,((double)N)/2.0))*exp(Blas_Dot_Prod(xMinusMean,invCovarianceXminusMean));
}

double StatUtil::Variance(const tVector &v)
{
	double squareMean=0.0,mean=0.0;

	for(int i=0;i<v.size();i++)
	{
		mean += v(i);
		squareMean += v(i)*v(i);
	}
	mean /= v.size();
	squareMean /= v.size();

	return squareMean - mean*mean;
}
