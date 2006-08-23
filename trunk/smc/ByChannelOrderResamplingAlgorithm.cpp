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
#include "ByChannelOrderResamplingAlgorithm.h"

ByChannelOrderResamplingAlgorithm::ByChannelOrderResamplingAlgorithm(ResamplingCriterion resamplingCriterion): ResamplingAlgorithm(resamplingCriterion)
{
}


ByChannelOrderResamplingAlgorithm::~ByChannelOrderResamplingAlgorithm()
{
}


void ByChannelOrderResamplingAlgorithm::Resample(ParticleFilter *particleFilter)
{
    ParticleFilterWithChannelOrder *pf = dynamic_cast <ParticleFilterWithChannelOrder *> (particleFilter);

    vector<vector <int> > indexesOfChannelOrders = pf->GetIndexesOfChannelOrders();
    tVector weights = pf->GetWeightsVector();

    for(int iChannelOrder=0;iChannelOrder<pf->NchannelOrders();iChannelOrder++)
    {
        int nParticlesCurrentChannelOrder = pf->NparticlesOfChannelOrderIndex(iChannelOrder);
        // if there is no particles with this channel order
        if(nParticlesCurrentChannelOrder==0)
            continue;

        if(_resamplingCriterion.ResamplingNeeded(weights,indexesOfChannelOrders[iChannelOrder]))
        {
            // the weights corresponding to the channel order being processed are selected. Note that they all are adjacent
//             tRange rWeightsOrder(indexesOfChannelOrders[iChannelOrder][0],indexesOfChannelOrders[iChannelOrder][indexesOfChannelOrders[iChannelOrder].size()-1]);

            tVector currentChannelOrderWeights(nParticlesCurrentChannelOrder);
            for(int i=0;i<nParticlesCurrentChannelOrder;i++)
                currentChannelOrderWeights(i) = weights(indexesOfChannelOrders[iChannelOrder][i]);
//
// //          cout << rWeightsOrder << endl;
//
            vector<int> auxIndexes = StatUtil::Discrete_rnd(nParticlesCurrentChannelOrder,currentChannelOrderWeights);
//
            vector<int> indexes(nParticlesCurrentChannelOrder);
            for(int i=0;i<nParticlesCurrentChannelOrder;i++)
                    indexes[i] = indexesOfChannelOrders[iChannelOrder][auxIndexes[i]];

            pf->KeepParticles(indexes,indexesOfChannelOrders[iChannelOrder]);
        }
    }
}

