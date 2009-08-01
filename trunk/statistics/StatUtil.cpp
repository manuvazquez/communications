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

#include <defines.h>

// #define DEBUG

using namespace std;

#ifdef RANDOM_SEED
Random StatUtil::_randomGenerator;
Random StatUtil::_particlesInitializerRandomGenerator;
#else
Random StatUtil::_randomGenerator(10);
Random StatUtil::_particlesInitializerRandomGenerator(20);
// Random StatUtil::_particlesInitializerRandomGenerator(200);
#endif

int StatUtil::discrete_rnd(const tVector &probabilities,Random &randomGenerator)
{
    int i;
    double uniform;

//     tVector normalizedProbabilities = Util::normalize(probabilities);
    int nProbabilities = probabilities.size();

    double *distributionFunction = new double[nProbabilities];
    distributionFunction[0] = probabilities(0);
    for(i=1;i<nProbabilities;i++)
           distributionFunction[i] = distributionFunction[i-1]+probabilities(i);

    uniform = randomGenerator.rand();
#ifdef DEBUG
    cout << "uniform es " << uniform << endl;
    cout << "distributionFunction[0] = " << distributionFunction[0] << endl;
#endif
    int res = 0;
    while(uniform>distributionFunction[res])
        res++;

    // memory release
    delete[] distributionFunction;

    return res;
}

int StatUtil::discrete_rnd(const std::vector<double> &probabilities,Random &randomGenerator)
{
    int i;
    double uniform;

    int nProbabilities = probabilities.size();

    double *distributionFunction = new double[nProbabilities];
    distributionFunction[0] = probabilities[0];
    for(i=1;i<nProbabilities;i++)
           distributionFunction[i] = distributionFunction[i-1]+probabilities[i];

    uniform = randomGenerator.rand();

    int res = 0;
    while(uniform>distributionFunction[res])
        res++;

    // memory release
    delete[] distributionFunction;

    return res;
}

vector<int> StatUtil::discrete_rnd(int nSamples,const tVector &probabilities,Random &randomGenerator)
{
    int i,j;
    double uniform;

    tVector normalizedProbabilities = Util::normalize(probabilities);
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

tMatrix StatUtil::randnMatrix(int rows,int cols,double mean,double variance,Random &randomGenerator)
{
    tMatrix res(rows,cols);
    double stdDv = sqrt(variance);

    int j;
    for(int i=0;i<rows;i++)
        for(j=0;j<cols;j++)
            res(i,j) =  randomGenerator.randn()*stdDv + mean;

    return res;
}

tVector StatUtil::randMatrix(const tVector &mean,const tMatrix &covariance,Random &randomGenerator)
{
    if(covariance.rows()!=mean.size() || covariance.cols()!=mean.size())
        throw RuntimeException("StatUtil::randnMatrix: dimensions of the mean or the covariance are wrong.");

    tVector res = mean;
    // res = mean + L*randnMatrix(mean.size(),1,0.0,1.0)
    Blas_Mat_Vec_Mult(Util::cholesky(covariance),randnMatrix(mean.size(),1,0.0,1.0,randomGenerator),res,1.0,1.0);

    return res;
}

double StatUtil::normalPdf(double x,double mean,double variance)
{
    double distance = fabs(x - mean);

    return 1.0/sqrt(2.0*M_PI*variance)*exp(-(distance*distance)/(2.0*variance));
}

double StatUtil::normalPdf(const tVector &x,const tVector &mean,const tMatrix &covariance)
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
    Util::add(x,mean,xMinusMean,1.0,-1.0);

    tVector invCovarianceXminusMean(N);
    // invCovarianceXminusMean = -0.5 * invCovariance * xMinusMean
    Blas_Mat_Vec_Mult(invCovariance,xMinusMean,invCovarianceXminusMean,-0.5);

    return 1.0/(sqrt(fabs(detCovariance))*pow(2.0*M_PI,((double)N)/2.0))*exp(Blas_Dot_Prod(xMinusMean,invCovarianceXminusMean));
}

double StatUtil::normalPdf(const tVector &x,const tVector &mean,double variance)
{
//  tMatrix covariance = LaGenMatDouble::eye(x.size(),x.size());
// 
//  covariance *= variance;
//  return StatUtil::normalPdf(x,mean,covariance);
    
    double res = 1.0;
    
    for(uint i=0;i<static_cast<uint> (x.rows());i++)
        res *= 1.0/sqrt(2.0*M_PI*variance)*exp(-((x(i) - mean(i))*(x(i) - mean(i)))/(2.0*variance));

    return res;
}

double StatUtil::variance(const tVector &v)
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

double StatUtil::mean(const tMatrix &A)
{
    double sum = 0.0;

    for(int i=0;i<A.rows();i++)
        for(int j=0;j<A.cols();j++)
            sum += A(i,j);

    return sum/(double)(A.rows()*A.cols());
}

vector<int> StatUtil::withoutReplacementSampling(int nSamples,const tVector &probabilities,Random &randomGenerator)
{
//     int i,j;
    double uniform;

    int nProbabilities = probabilities.size();

    if(nSamples>=nProbabilities)
    {
        vector<int> res(nProbabilities);
        for(int i=0;i<nProbabilities;i++)
            res[i] = i;
        return res;
    }

    tVector normalizedProbabilities = Util::normalize(probabilities);

    bool **distributionFunctionActiveOperands;
    distributionFunctionActiveOperands = new bool*[nProbabilities];

    bool *remainingProbabilityActiveOperands;
    remainingProbabilityActiveOperands = new bool[nProbabilities];
    for(int i=0;i<nProbabilities;i++)
    {
        distributionFunctionActiveOperands[i] = new bool[nProbabilities];
        for(int j=0;j<=i;j++)
        {
            distributionFunctionActiveOperands[i][j] = true;
            remainingProbabilityActiveOperands[j] = true;
        }
        for(int j=i+1;j<nProbabilities;j++)
        {
            distributionFunctionActiveOperands[i][j] = false;
            remainingProbabilityActiveOperands[j] = true;
        }
    }

    int j;
    vector<int> res(nSamples);
    for(int i=0;i<nSamples;i++)
    {
        uniform = randomGenerator.rand()*computeFromActiveOperands(normalizedProbabilities,remainingProbabilityActiveOperands);

        j=0;
        while(uniform > computeFromActiveOperands(normalizedProbabilities,distributionFunctionActiveOperands[j]))
            j++;
        res[i] = j;

        remainingProbabilityActiveOperands[j] = false;
        for(int k=j;k<nProbabilities;k++)
            distributionFunctionActiveOperands[k][j] = false;
    }

    for(int i=0;i<nProbabilities;i++)
        delete[] distributionFunctionActiveOperands[i];
    delete[] distributionFunctionActiveOperands;
    delete[] remainingProbabilityActiveOperands;

    return res;
}

inline double StatUtil::computeFromActiveOperands(const tVector &probabilities,bool *activeOperands)
{
    double res = 0.0;
    for(int i=0;i<probabilities.size();i++)
        if(activeOperands[i])
            res += probabilities(i);
    return res;
}
