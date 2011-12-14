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

// the seed used to create the random objects is generated from the system time
//  #define RANDOM_SEED

// whether the ENTIRE sequence of channel matrices estimated by every path is kept (it doesn't get saved, either way)
// this only applies to algorithms using "PSPPath"
// #define DO_NOT_STORE_THE_SEQUENCE_OF_CHANNEL_MATRICES_ESTIMATED_BY_EVERY_PATH

#define ALGORITHM_NAME_MAX_LENGTH 80

// used as a value when there shouldn't be any value :) (e.g. not computed BER's or MSE's)
#define FUNNY_VALUE -3.14
