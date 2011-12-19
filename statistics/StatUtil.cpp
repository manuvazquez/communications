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

// #ifdef RANDOM_SEED
  Random StatUtil::_randomGenerator;
  Random StatUtil::_particlesInitializerRandomGenerator;
// #else
//   Random StatUtil::_randomGenerator(4135925433);
//   Random StatUtil::_particlesInitializerRandomGenerator(2484546298);
// #endif

uint StatUtil::discrete_rnd(const VectorXd &probabilities,Random &randomGenerator)
{
    uint i;
    double uniform;

    uint nProbabilities = probabilities.size();

    double *distributionFunction = new double[nProbabilities];
    distributionFunction[0] = probabilities(0);
    for(i=1;i<nProbabilities;i++)
           distributionFunction[i] = distributionFunction[i-1]+probabilities(i);

    uniform = randomGenerator.rand();
    
    uint res = 0;
    while(uniform>distributionFunction[res])
        res++;

    // memory release
    delete[] distributionFunction;

    return res;
}

uint StatUtil::discrete_rnd(const std::vector<double> &probabilities,Random &randomGenerator)
{
    double uniform;

    uint nProbabilities = probabilities.size();

    double *distributionFunction = new double[nProbabilities];
    distributionFunction[0] = probabilities[0];
    for(uint i=1;i<nProbabilities;i++)
           distributionFunction[i] = distributionFunction[i-1]+probabilities[i];

    uniform = randomGenerator.rand();

    uint res = 0;
    while(uniform>distributionFunction[res])
        res++;

    // memory release
    delete[] distributionFunction;

    return res;
}

vector<uint> StatUtil::discrete_rnd(uint nSamples,const VectorXd &probabilities,Random &randomGenerator)
{
    uint i,j;
    double uniform;

    VectorXd normalizedProbabilities = Util::normalize(probabilities);
    uint nProbabilities = probabilities.size();

    double *distributionFunction = new double[nProbabilities];
    distributionFunction[0] = normalizedProbabilities(0);
    for(i=1;i<nProbabilities;i++)
           distributionFunction[i] = distributionFunction[i-1]+normalizedProbabilities(i);

    vector<uint> res(nSamples);

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

MatrixXd StatUtil::randnMatrix(uint rows,uint cols,double mean,double variance,Random &randomGenerator)
{
    MatrixXd res(rows,cols);
    double stdDv = sqrt(variance);

    uint j;
    for(uint i=0;i<rows;i++)
        for(j=0;j<cols;j++)
            res(i,j) =  randomGenerator.randn()*stdDv + mean;

    return res;
}

VectorXd StatUtil::randnMatrix(const VectorXd &mean,const MatrixXd &covariance,Random &randomGenerator)
{
    if(covariance.rows()!=mean.size() || covariance.cols()!=mean.size())
        throw RuntimeException("StatUtil::randnMatrix: dimensions of the mean or the covariance are wrong.");

//     return mean + Util::cholesky(covariance)*randnMatrix_eigen(mean.size(),1,0.0,1.0,randomGenerator);
    return mean + Eigen::LLT<MatrixXd>(covariance).matrixL()*randnMatrix(mean.size(),1,0.0,1.0,randomGenerator);
}

double StatUtil::normalPdf(double x,double mean,double variance)
{
    double distance = fabs(x - mean);

    return 1.0/sqrt(2.0*M_PI*variance)*exp(-(distance*distance)/(2.0*variance));
}

double StatUtil::normalPdf(const VectorXd &x,const VectorXd &mean,const MatrixXd &covariance)
{
    uint N = x.size();
    
    Eigen::LDLT<MatrixXd> ldltOfCovariance(covariance);

    MatrixXd invCovariance = MatrixXd::Identity(N,N);
    ldltOfCovariance.solveInPlace(invCovariance);
        
    double invCovarianceDeterminant = 1.0;    
    for(uint i=0;i<ldltOfCovariance.vectorD().rows();i++)
        invCovarianceDeterminant *= ldltOfCovariance.vectorD().coeff(i); 

    return 1.0/(sqrt(fabs(invCovarianceDeterminant))*pow(2.0*M_PI,((double)N)/2.0))*exp((x-mean).dot(-0.5*invCovariance*(x-mean)));
}

double StatUtil::normalPdf(const VectorXd &x,const VectorXd &mean,double variance)
{
    double res = 1.0;
    
    for(uint i=0;i<x.rows();i++)
        res *= 1.0/sqrt(2.0*M_PI*variance)*exp(-((x(i) - mean(i))*(x(i) - mean(i)))/(2.0*variance));

    return res;
}

double StatUtil::variance(const VectorXd &v)
{
    double squareMean=0.0,mean=0.0;

    for(uint i=0;i<v.size();i++)
    {
        mean += v(i);
        squareMean += v(i)*v(i);
    }
    mean /= v.size();
    squareMean /= v.size();

    return squareMean - mean*mean;
}

double StatUtil::mean(const MatrixXd &A)
{
    double sum = 0.0;

    for(uint i=0;i<A.rows();i++)
        for(uint j=0;j<A.cols();j++)
            sum += A(i,j);

    return sum/(double)(A.rows()*A.cols());
}

vector<uint> StatUtil::withoutReplacementSampling(uint nSamples, const VectorXd& probabilities, Random& randomGenerator)
{
    double uniform;

    uint nProbabilities = probabilities.size();

    if(nSamples>=nProbabilities)
    {
        vector<uint> res(nProbabilities);
        for(uint i=0;i<nProbabilities;i++)
            res[i] = i;
        return res;
    }

    VectorXd normalizedProbabilities = Util::normalize(probabilities);

    bool **distributionFunctionActiveOperands;
    distributionFunctionActiveOperands = new bool*[nProbabilities];

    bool *remainingProbabilityActiveOperands;
    remainingProbabilityActiveOperands = new bool[nProbabilities];
    for(uint i=0;i<nProbabilities;i++)
    {
        distributionFunctionActiveOperands[i] = new bool[nProbabilities];
        for(uint j=0;j<=i;j++)
        {
            distributionFunctionActiveOperands[i][j] = true;
            remainingProbabilityActiveOperands[j] = true;
        }
        for(uint j=i+1;j<nProbabilities;j++)
        {
            distributionFunctionActiveOperands[i][j] = false;
            remainingProbabilityActiveOperands[j] = true;
        }
    }

    uint j;
    vector<uint> res(nSamples);
    for(uint i=0;i<nSamples;i++)
    {
        uniform = randomGenerator.rand()*computeFromActiveOperands(normalizedProbabilities,remainingProbabilityActiveOperands);

        j=0;
        while(uniform > computeFromActiveOperands(normalizedProbabilities,distributionFunctionActiveOperands[j]))
            j++;
        res[i] = j;

        remainingProbabilityActiveOperands[j] = false;
        for(uint k=j;k<nProbabilities;k++)
            distributionFunctionActiveOperands[k][j] = false;
    }

    for(uint i=0;i<nProbabilities;i++)
        delete[] distributionFunctionActiveOperands[i];
    delete[] distributionFunctionActiveOperands;
    delete[] remainingProbabilityActiveOperands;

    return res;
}

inline double StatUtil::computeFromActiveOperands(const VectorXd &probabilities,bool *activeOperands)
{
    double res = 0.0;
    for(uint i=0;i<probabilities.size();i++)
        if(activeOperands[i])
            res += probabilities(i);
    return res;
}

double StatUtil::probApriori(const VectorXd &symbolsVector, const std::vector<UsersActivityDistribution> &symbolsDistributions)
{
  if(static_cast<uint>(symbolsVector.size())!=symbolsDistributions.size())
	throw RuntimeException("StatUtil::probApriori: the number of symbols in the vector and that of distributions don't match.");
  
  double res = 1.0;
  
  for(uint i=0;i<symbolsVector.size();i++)
	res *= symbolsDistributions[i].probApriori(Util::isUserActive(symbolsVector(i)));
// 	res *= symbolsDistributions[i].probApriori(symbolsDistributions[i].isUserActive(symbolsVector(i)));  
  
  return res;
}

double StatUtil::probXgivenY(VectorXd &X, VectorXd &Y, const std::vector<UsersActivityDistribution> &symbolsDistributions)
{ 
  if(X.size()!=Y.size())
	throw RuntimeException("StatUtil::probXgivenY: the sizes of the vectors don't match.");

  if(static_cast<uint>(X.size())!=symbolsDistributions.size())
	throw RuntimeException("StatUtil::probXgivenY: the number of symbols in the vectors and that of distributions don't match.");  
  
  double res = 1.0;
  
  for(uint i=0;i<X.size();i++)
	res *= symbolsDistributions[i].probXgivenY(Util::isUserActive(X(i)),Util::isUserActive(Y(i)));
// 	res *= symbolsDistributions[i].probXgivenY(symbolsDistributions[i].isUserActive(X(i)),symbolsDistributions[i].isUserActive(Y(i)));  
  
  return res;
}

// this definitely needs some optimizing
double StatUtil::probSymbolsVectorGivenPreviousTimeInstantUsersActivity(const VectorXd& symbolsVector, const std::vector< bool >& previousTimeInstantUsersActivity, const std::vector<UsersActivityDistribution> &usersActivityPdfs, uint alphabetLength)
{
    if(static_cast<uint> (symbolsVector.size())!=previousTimeInstantUsersActivity.size())
        throw RuntimeException("StatUtil::probSymbolsVectorGivenPreviousTimeInstantUsersActivity: symbols vector size doesn't coincide with that of the vector containing information about the users activity in the previous time instant.");
        
    double probSymbolWhenUserActive,probSymbolWhenUserNotActive;
	double overallProb = 1.0;
    
    for(uint i=0;i<symbolsVector.size();i++)
    {
	  if(Util::isUserActive(symbolsVector(i)))
	  {
		probSymbolWhenUserNotActive = 0.0;
		probSymbolWhenUserActive = 1/double(alphabetLength) * usersActivityPdfs[i].probXgivenY(true,previousTimeInstantUsersActivity[i]);
	  }
	  else
	  {
		probSymbolWhenUserNotActive = usersActivityPdfs[i].probXgivenY(false,previousTimeInstantUsersActivity[i]);
		probSymbolWhenUserActive = 0.0;
	  }

	  overallProb *= (probSymbolWhenUserNotActive + probSymbolWhenUserActive);
    }

    return overallProb;
}

double StatUtil::probSymbolsVector(const VectorXd &symbolsVector,const std::vector<UsersActivityDistribution> &usersActivityPdfs, uint alphabetLength)
{
  double probSymbolWhenUserActive,probSymbolWhenUserNotActive;
  double overallProb = 1.0;

  for(uint i=0;i<symbolsVector.size();i++)
  {
	if(Util::isUserActive(symbolsVector(i)))
	{
	  probSymbolWhenUserNotActive = 0.0;
	  probSymbolWhenUserActive = 1/double(alphabetLength) * usersActivityPdfs[i].probApriori(true);
	}
	else
	{
	  probSymbolWhenUserNotActive = usersActivityPdfs[i].probApriori(false);
	  probSymbolWhenUserActive = 0.0;
	}

	overallProb *= (probSymbolWhenUserNotActive + probSymbolWhenUserActive);
  }

  return overallProb;
}