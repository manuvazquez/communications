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
#ifndef PARTICLE_H
#define PARTICLE_H

/**
    @author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <Util.h>
#include <types.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>

class Particle{
protected:
    double _weight;
    MatrixXd _symbolVectors;
public:
    Particle(double weight,int symbolVectorLength,int nTimeInstants);
    virtual ~Particle();

    /**
     *
     * @return the number of time instants
     */
    int trajectorylength() const { return _symbolVectors.cols();}

    double getWeight() const { return _weight;}
    void setWeight(double weight) { _weight = weight;}

    tMatrix getAllSymbolVectors() const { return Util::eigen2lapack(_symbolVectors);}
    
    tVector getSymbolVector(int n) const
    { 
        VectorXd aux = _symbolVectors.col(n);
        return Util::eigen2lapack(aux);
    }
    VectorXd getSymbolVector_eigen(int n) const { return _symbolVectors.col(n);}
    
    void setSymbolVector(int n,const tVector &v) { _symbolVectors.col(n) = Util::lapack2eigen(v);}
    void setSymbolVector(int n,const VectorXd &v) { _symbolVectors.col(n) = v;}
    void setSymbolVector(int n,const std::vector<tSymbol> &v)
    {
        for(int i=0;i<_symbolVectors.rows();i++)
            _symbolVectors(i,n) = v[i];
    }

    tMatrix getSymbolVectors(const tRange &range) const { return Util::eigen2lapack(_symbolVectors)(tRange(0,_symbolVectors.rows()-1),range);}
    tMatrix getSymbolVectors(int a,int b) const { return Util::eigen2lapack(_symbolVectors)(tRange(0,_symbolVectors.rows()-1),tRange(a,b));}
    MatrixXd getSymbolVectors() { return _symbolVectors;}

    void setSymbolVectors(const tRange &range,const tMatrix &symbolVectors)
    {
        setSymbolVectors(range.start(),range.end()+1,Util::lapack2eigen(symbolVectors));
    }

    void setSymbolVectors(int a,int b,const tMatrix &symbolVectors)
    {
        setSymbolVectors(a,b+1,Util::lapack2eigen(symbolVectors));
    }

    void setSymbolVectors(int a,int b,const MatrixXd &symbolVectors)
    {
        _symbolVectors.block(0,a,_symbolVectors.rows(),b-a) = symbolVectors;
    }

    void print() const { std::cout << _symbolVectors << std::endl << "peso = " << _weight << std::endl;}

    virtual Particle *clone()
    {
        return new Particle(*this);
    }

    void operator=(const Particle &particle);
};

#endif
