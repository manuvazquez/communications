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
#include <TesisComplejidadReducidaSystem.h>
#include <TesisComplejidadReducidaBesselSystem.h>
#include <TesisComplejidadReducidaARSystem.h>
#include <TesisOrdenCanalDesconocidoARSystem.h>
#include <TesisOrdenCanalDesconocidoBesselSystem.h>
#include <TesisComplejidadReducidaBesselNumeroParticulasSystem.h>
#include <LMSmuTestSystem.h>
#include <PSPvsSMCSystem.h>
#include <CDMASystem.h>
#include <ISWCS10System.h>

#include <signal.h>

bool __done = false;
bool __randomSeedHasBeenPassed = false;
bool __nFramesHasBeenPassed = false;

uint __nFramesPassed;
uint32_t __mainSeedPassed;
uint32_t __statUtilSeedPassed;

void setDoneTrue(int signal)
{
	std::cout << "Ctl+C read. Finishing frame..." << std::endl;
	__done  = true;
}

int main(int argc,char* argv[])
{
// 	signal(SIGINT,&setDoneTrue);
// 
// 	Elsevier2007BesselChannelSystem system;
// 	Elsevier2007ARChannelSystem  system;
// 	TVT2007System system;
// 	WSA08System system;
// 	PSPvsPSPBasedSMCSystem system;
// 	Rev2TVT2007System system;
// 
// 	LMSmuTestSystem system;
// 	PSPvsSMCSystem system;
// 
// 	TesisOrdenCanalDesconocidoARSystem system;
// 	TesisOrdenCanalDesconocidoBesselSystem system;
// 
// 	TesisComplejidadReducidaBesselSystem system;
// 	TesisComplejidadReducidaARSystem system;
// 
// 	TesisComplejidadReducidaBesselNumeroParticulasSystem system;
// 
// 	CDMASystem system;

	std::cout << "received " << argc << " arguments" << endl;
	
// 	if(argc>1)
// 	{
// 		// seed (char *) is converted to a number
// 		std::istringstream stream(argv[1]);
// 		
// 		// to check if everything went well: if (stream.eof())
// 
// 		stream >> __randomSeedPassed;
// 		cout << "received seed is " << __randomSeedPassed << endl;
// 		__randomSeedHasBeenPassed = true;
// 	}else 

	std::istringstream argument;

	switch(argc)
	{
		case 1:
			cout << "no arguments received" << endl;
			break;
			
		case 2:
			argument.clear();
			argument.str(argv[1]);
			
			// to check if everything went well: if (argument.eof())

			argument >> __mainSeedPassed;
			
			// we assume both seeds are identical
			__statUtilSeedPassed = __mainSeedPassed;
			
			__randomSeedHasBeenPassed = true;
			cout << "received seed is " << __mainSeedPassed << endl;
			break;
			
		case 4:
			argument.clear();
			argument.str(argv[1]);
			argument >> __mainSeedPassed;
			
			argument.clear();
			argument.str(argv[2]);
			argument >> __statUtilSeedPassed;
			
			argument.clear();
			argument.str(argv[3]);
			argument >> __nFramesPassed;
			
			__randomSeedHasBeenPassed = true;
			__nFramesHasBeenPassed = true;
			cout << "assuming:" << endl << "\t" << "mainSeed = " << __mainSeedPassed << endl << "\t" << "statUtilSeed = " << __statUtilSeedPassed << endl<< "\t" << "nFrames = " << __nFramesPassed << endl;
			break;
			
		default:
			cout << "incorrect number of arguments!!" << endl;
			exit(1);
	}
		

	ISWCS10System system;

	system.simulate();
}
