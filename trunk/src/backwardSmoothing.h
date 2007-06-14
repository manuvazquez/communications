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
vector<int> testingBackwardsSmoothings;

testingBackwardsSmoothings.push_back(5);
testingBackwardsSmoothings.push_back(3);
testingBackwardsSmoothings.push_back(0);
vector<LinearDetector *> RMMSEbackwardsSmoothingTesting;

for(uint iBackwardsSmoothing=0;iBackwardsSmoothing<testingBackwardsSmoothings.size();iBackwardsSmoothing++)
	RMMSEbackwardsSmoothingTesting.push_back(new RMMSEDetector(L*(testingBackwardsSmoothings[iBackwardsSmoothing]+m),N*(m+testingBackwardsSmoothings[iBackwardsSmoothing]+m-1),pam2.Variance(),forgettingFactorDetector,N*m));
#else
// the RLS-D-SIS algorithm for serveral backward smoothing parameters
for(uint iBackwardsSmoothing=0;iBackwardsSmoothing<testingBackwardsSmoothings.size();iBackwardsSmoothing++)
{
	char buffer[SPRINTF_BUFFER];

	sprintf(buffer," backwards smoothing = %d",testingBackwardsSmoothings[iBackwardsSmoothing]);

	algorithms.push_back(new LinearFilterBasedMKFAlgorithm(string("MKF") + string(buffer),pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,RMMSEbackwardsSmoothingTesting[iBackwardsSmoothing],preamble,testingBackwardsSmoothings[iBackwardsSmoothing],d,nParticles,&algoritmoRemuestreo,initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
}
#endif
