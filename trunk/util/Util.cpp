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
#include "Util.h"

void Util::Add(const tMatrix& A,const tMatrix& B,tMatrix& C,double alpha,double beta)
{
	int i,j;
	int rows = A.rows(), cols = A.cols();
	for(i=0;i<rows;i++)
		for(j=0;j<cols;j++)
			C(i,j) = alpha*A(i,j) + beta*B(i,j);
}

// void Util::GaussianMatrix(tMatrix& matrix,double mean,double variance)
// {
//
// }
