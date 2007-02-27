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
public:
	double _cost;
	tMatrix *_detectedSequence;
public:
    ViterbiPath();
    ViterbiPath(const ViterbiPath &path);
    ViterbiPath(const ViterbiPath &path, tVector newSymbolVector, double newCost);
    ~ViterbiPath();

	void Update(const ViterbiPath &path, tVector newSymbolVector, double newCost);
	void Ref(const ViterbiPath &path);
	void operator=(const ViterbiPath &path);
};

#endif
