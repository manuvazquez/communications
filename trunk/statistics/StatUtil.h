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
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>
#include <lapackpp/sybmd.h>

#include <Eigen/Cholesky>

class StatUtil{
private:
    static double computeFromActiveOperands(const tVector &probabilities,bool *activeOperands);
    static Random _randomGenerator;
public:
    static Random _particlesInitializerRandomGenerator;


    /**
     * It assumes that the probabilities are normalized
     * @param probabilities
     * @return
     */
    static int discrete_rnd(const tVector &probabilities,Random &randomGenerator = _randomGenerator);
    static int discrete_rnd(const VectorXd &probabilities,Random &randomGenerator = _randomGenerator); // eigen
    static vector<int> discrete_rnd(int nSamples,const tVector &probabilities,Random &randomGenerator = _randomGenerator);
    
    // same functions as above but receiving c++ vectors instead of Lapack++ vectors
    static int discrete_rnd(const std::vector<double> &probabilities,Random &randomGenerator = _randomGenerator);    
    
    static tMatrix randnMatrix(int rows,int cols,double mean,double variance,Random &randomGenerator = _randomGenerator);
    static MatrixXd randnMatrix_eigen(int rows,int cols,double mean,double variance,Random &randomGenerator = _randomGenerator);
    static tVector randMatrix(const tVector &mean,const tMatrix &covariance,Random &randomGenerator = _randomGenerator);
    static VectorXd randnMatrix(const VectorXd &mean,const MatrixXd &covariance,Random &randomGenerator = _randomGenerator); // eigen
    static double normalPdf(double x,double mean,double variance);
    static double normalPdf(const tVector &x,const tVector &mean,const tMatrix &covariance);
    static double normalPdf(const tVector &x,const tVector &mean,double variance);
    static double normalPdf(const VectorXd &x,const VectorXd &mean,const MatrixXd &covariance); // eigen
    static double normalPdf(const VectorXd &x,const VectorXd &mean,double variance); //eigen
    static double variance(const tVector &v);
    static double mean(const tMatrix &A);
    static vector<int> withoutReplacementSampling(int nSamples,const tVector &probabilities,Random &randomGenerator = _randomGenerator);
    static Random& getRandomGenerator() { return _randomGenerator;}
};

#endif
