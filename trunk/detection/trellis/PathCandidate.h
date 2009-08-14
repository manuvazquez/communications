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
#ifndef PATHCANDIDATE_H
#define PATHCANDIDATE_H

/**
	@author Manu <manu@rustneversleeps>
*/
class PathCandidate{
public:
	int _fromState;
	int _input;
	double _cost;
// 	tVector _newSymbolVector;
    VectorXd _newSymbolVector;

    PathCandidate():_cost(-1.0) {}

    void Clean() { _cost = -1.0;}
    double GetCost() { return _cost;}
    bool NoPathArrived() { return (_cost < 0);}
};

#endif
