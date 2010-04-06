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
#ifndef VITERBIPATH_H
#define VITERBIPATH_H

/**
    @author Manu <manu@rustneversleeps>
*/

#include <types.h>
#include <exceptions.h>

class ViterbiPath{
protected:
    int _nTimeInstants;
    double _cost;
    MatrixXd *_detectedSequence;
public:
    ViterbiPath();
    ViterbiPath(int nTimeInstants,double cost,MatrixXd initialSequence);
	ViterbiPath(int nTimeInstants,double cost);
    ViterbiPath(const ViterbiPath &path);
    virtual ~ViterbiPath();

    double getCost() const { return _cost;}
    void clean() { delete _detectedSequence; _detectedSequence = NULL;}
    bool isEmpty() const { return (_detectedSequence == NULL);}
    VectorXd getSymbolVector(int n) const { return _detectedSequence->col(n);}
    MatrixXd getDetectedSequence()  const { return *_detectedSequence;}

    virtual void print() const;
	
	/*! 
	  It updates this path by overwriting it with a new one built from \ref path
	  \param path the other \ref ViterbiPath from which most of the information to build the new one will be taken
	  \param newSymbolVector this will be added to the sequence accumulated by \ref path
	  \param newCost the cost assigned to the new \ref ViterbiPath
	*/
    void update(const ViterbiPath &path, VectorXd newSymbolVector, double newCost); // eigen
    
	virtual void Ref(const ViterbiPath &path);
    virtual void operator=(const ViterbiPath &path);
};

#endif
