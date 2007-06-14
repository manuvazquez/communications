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
vector<double> resamplingRates;
resamplingRates.push_back(0.001);
resamplingRates.push_back(0.05);
resamplingRates.push_back(0.1);
resamplingRates.push_back(0.3);
resamplingRates.push_back(0.5);
resamplingRates.push_back(0.9);

vector<ResidualResamplingAlgorithm> resamplingAlgorithms;
for(uint iResamplingAlgorithm=0;iResamplingAlgorithm<resamplingRates.size();iResamplingAlgorithm++)
	resamplingAlgorithms.push_back(ResidualResamplingAlgorithm(ResamplingCriterion(resamplingRates[iResamplingAlgorithm])));
#else
// the RLS algorithm considering differente resampling rates
for(uint iResamplingAlgorithm=0;iResamplingAlgorithm<resamplingRates.size();iResamplingAlgorithm++)
{
	char buffer[SPRINTF_BUFFER];

	sprintf(buffer," resampling rate = %f",resamplingRates[iResamplingAlgorithm]);

	algorithms.push_back(new LinearFilterBasedSMCAlgorithm(string("RLS-D-SIS") + string(buffer),pam2,L,N,lastSymbolVectorInstant,m,&rlsEstimator,&rmmseDetector,preamble,0,d,nParticles,&(resamplingAlgorithms[iResamplingAlgorithm]),initialChannelEstimation,channelCoefficientsVariances,ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
}
#endif
