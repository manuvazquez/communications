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

#ifndef PARAMETERS_DEFINED
vector<int> __particlesNumberTesting;
__particlesNumberTesting.push_back(30);__particlesNumberTesting.push_back(50);
__particlesNumberTesting.push_back(100);__particlesNumberTesting.push_back(200);
__particlesNumberTesting.push_back(500);__particlesNumberTesting.push_back(1000);
#else
for(uint i=0;i<__particlesNumberTesting.size();i++)
{
	char buffer[SPRINTF_BUFFER];

	sprintf(buffer," (%d particles)",__particlesNumberTesting[i]);

	algorithms.push_back(new LinearFilterBasedSMCAlgorithm(string("RLS-D-SIS") + string(buffer),pam2,L,N,lastSymbolVectorInstant,m,&rlsEstimator,&rmmseDetector,preamble,0,d,__particlesNumberTesting[i],&algoritmoRemuestreo,powerProfile.Means(),powerProfile.Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
}
#endif
