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
#ifndef STATUTIL_H
#define STATUTIL_H

/**
@author Manu
*/

#include <stdlib.h>
#include <vector>
#include <types.h>
#include <Random.h>
#include <Util.h>
#include <UsersActivityDistribution.h>

#include <Eigen/Cholesky>

class StatUtil{
private:
    static double computeFromActiveOperands(const VectorXd &probabilities,bool *activeOperands); // eigen
    static Random _randomGenerator;
public:
	// this is only used in USIS =>
	// FIXME: it should be taken out
    static Random _particlesInitializerRandomGenerator;


    /**
     * It assumes that the probabilities are normalized
     * @param probabilities
     * @return
     */
    static uint discrete_rnd(const VectorXd &probabilities,Random &randomGenerator = _randomGenerator); // eigen
    static vector< uint > discrete_rnd(uint nSamples, const VectorXd& probabilities, Random& randomGenerator = _randomGenerator); // eigen
    
    // same functions as above but receiving c++ vectors instead of "eigen" vectors
    static uint discrete_rnd(const std::vector<double> &probabilities,Random &randomGenerator = _randomGenerator);    
    
    static MatrixXd randnMatrix(int rows,int cols,double mean,double variance,Random &randomGenerator = _randomGenerator);
    static VectorXd randnMatrix(const VectorXd &mean,const MatrixXd &covariance,Random &randomGenerator = _randomGenerator); // eigen
    static double normalPdf(double x,double mean,double variance);
    static double normalPdf(const VectorXd &x,const VectorXd &mean,const MatrixXd &covariance); // eigen
    static double normalPdf(const VectorXd &x,const VectorXd &mean,double variance); //eigen
    static double variance(const VectorXd &v);
    static double mean(const MatrixXd &A);
    static vector<uint> withoutReplacementSampling(uint nSamples,const VectorXd &probabilities,Random &randomGenerator = _randomGenerator); // eigen
    static Random& getRandomGenerator() { return _randomGenerator;}
    static double probApriori(const VectorXd &symbolsVector, const std::vector<UsersActivityDistribution> &symbolsDistributions);
    static double probXgivenY(VectorXd &X, VectorXd &Y, const std::vector<UsersActivityDistribution> &symbolsDistributions);
	
	static double probSymbolsVectorGivenPreviousTimeInstantUsersActivity(const VectorXd& symbolsVector, const std::vector< bool >& previousTimeInstantUsersActivity, const std::vector<UsersActivityDistribution> &usersActivityPdfs, uint alphabetLength);
	static double probSymbolsVector(const VectorXd &symbolsVector,const std::vector<UsersActivityDistribution> &usersActivityPdfs, uint alphabetLength);	
};

#endif
