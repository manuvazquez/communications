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
vector<int> __backwardSmoothing;
vector<int> __forwardSmoothing;

__backwardSmoothing.push_back(5);__forwardSmoothing.push_back(5);
__backwardSmoothing.push_back(3);__forwardSmoothing.push_back(3);
__backwardSmoothing.push_back(2);__forwardSmoothing.push_back(2);
__backwardSmoothing.push_back(0);__forwardSmoothing.push_back(1);
__backwardSmoothing.push_back(0);__forwardSmoothing.push_back(2);
__backwardSmoothing.push_back(0);__forwardSmoothing.push_back(3);
vector<LinearDetector *> backwardForwardSmoothingDetectors;

for(uint iSmoothing=0;iSmoothing<__backwardSmoothing.size();iSmoothing++)
	backwardForwardSmoothingDetectors.push_back(new MMSEDetector(L*(__backwardSmoothing[iSmoothing]+__forwardSmoothing[iSmoothing]+1),N*(m+__backwardSmoothing[iSmoothing]+__forwardSmoothing[iSmoothing]),pam2.Variance(),N*(__forwardSmoothing[iSmoothing]+1)));
// 	backwardForwardSmoothingDetectors.push_back(new RMMSEDetector(L*(__backwardSmoothing[iSmoothing]+__forwardSmoothing[iSmoothing]+1),N*(m+__backwardSmoothing[iSmoothing]+__forwardSmoothing[iSmoothing]),pam2.Variance(),forgettingFactorDetector,N*(__forwardSmoothing[iSmoothing]+1)));
#else
// the RLS-D-SIS algorithm for serveral backward smoothing parameters
for(uint iSmoothing=0;iSmoothing<__backwardSmoothing.size();iSmoothing++)
{
	char buffer[SPRINTF_BUFFER];

	sprintf(buffer," c = %d, d = %d",__backwardSmoothing[iSmoothing],__forwardSmoothing[iSmoothing]);

	algorithms.push_back(new LinearFilterBasedSMCAlgorithm(string("Linear Filter Based (Known Channel)") + string(buffer),pam2,L,N,lastSymbolVectorInstant,m,&knownChannelEstimator,backwardForwardSmoothingDetectors[iSmoothing],preamble,__backwardSmoothing[iSmoothing],__forwardSmoothing[iSmoothing],nParticles,&algoritmoRemuestreo,powerProfile.Means(),channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));

// 	algorithms.push_back(new LinearFilterBasedMKFAlgorithm(string("MKF") + string(buffer),pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,backwardForwardSmoothingDetectors[iSmoothing],preamble,__backwardSmoothing[iSmoothing],__forwardSmoothing[iSmoothing],nParticles,&algoritmoRemuestreo,powerProfile.Means(),channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
}
#endif
