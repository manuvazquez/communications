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

#include <SMCSystem.h>
#include <Elsevier2007BesselChannelSystem.h>
#include <Elsevier2007ARChannelSystem.h>
#include <TVT2007System.h>
#include <PSPvsPSPBasedSMCSystem.h>
#include <WSA08System.h>
#include <Rev2TVT2007System.h>
#include <TesisOrdenCanalSystem.h>
#include <TesisOrdenCanalMedianteSISSystem.h>
#include <TesisComplejidadReducidaSystem.h>
#include <TesisComplejidadReducidaBesselSystem.h>
#include <TesisComplejidadReducidaARSystem.h>
#include <TesisOrdenCanalDesconocidoARSystem.h>
#include <TesisOrdenCanalDesconocidoBesselSystem.h>
#include <LMSmuTestSystem.h>
#include <PSPvsSMCSystem.h>

#include <signal.h>

bool __done = false;

void setDoneTrue(int signal)
{
	std::cout << "Ctl+C read. Finishing frame..." << std::endl;
	__done  = true;
}

int main(int argc,char* argv[])
{
// 	signal(SIGINT,&setDoneTrue);

//     Elsevier2007BesselChannelSystem system;
//     Elsevier2007ARChannelSystem  system;
//     TVT2007System system;
//     WSA08System system;
//     PSPvsPSPBasedSMCSystem system;
//     Rev2TVT2007System system;

//     TesisOrdenCanalSystem system;
//     TesisOrdenCanalMedianteSISSystem system;

//     TesisOrdenCanalDesconocidoARSystem system;
    TesisOrdenCanalDesconocidoBesselSystem system;

//     TesisComplejidadReducidaBesselSystem system;
//     TesisComplejidadReducidaARSystem system;

//     LMSmuTestSystem system;
//     PSPvsSMCSystem system;
    system.Simulate();
}
