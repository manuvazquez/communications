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
vector<int> backwardSmoothing;
vector<int> forwardSmoothing;

backwardSmoothing.push_back(5);forwardSmoothing.push_back(5);
backwardSmoothing.push_back(3);forwardSmoothing.push_back(3);
backwardSmoothing.push_back(2);forwardSmoothing.push_back(2);
backwardSmoothing.push_back(0);forwardSmoothing.push_back(1);
backwardSmoothing.push_back(0);forwardSmoothing.push_back(2);
backwardSmoothing.push_back(0);forwardSmoothing.push_back(3);
vector<LinearDetector *> RMMSEbackwardForwardSmoothingTest;

for(uint iSmoothing=0;iSmoothing<backwardSmoothing.size();iSmoothing++)
	RMMSEbackwardForwardSmoothingTest.push_back(new RMMSEDetector(L*(backwardSmoothing[iSmoothing]+forwardSmoothing[iSmoothing]+1),N*(m+backwardSmoothing[iSmoothing]+forwardSmoothing[iSmoothing]),pam2.Variance(),forgettingFactorDetector,N*(forwardSmoothing[iSmoothing]+1)));
#else
// the RLS-D-SIS algorithm for serveral backward smoothing parameters
for(uint iSmoothing=0;iSmoothing<backwardSmoothing.size();iSmoothing++)
{
	char buffer[SPRINTF_BUFFER];

	sprintf(buffer," c = %d, d = %d",backwardSmoothing[iSmoothing],forwardSmoothing[iSmoothing]);

	algorithms.push_back(new LinearFilterBasedMKFAlgorithm(string("MKF") + string(buffer),pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,RMMSEbackwardForwardSmoothingTest[iSmoothing],preamble,backwardSmoothing[iSmoothing],forwardSmoothing[iSmoothing],nParticles,&algoritmoRemuestreo,initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
}
#endif
