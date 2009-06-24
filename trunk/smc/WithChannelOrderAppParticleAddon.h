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
#ifndef WITHCHANNELORDERAPPPARTICLEADDON_H
#define WITHCHANNELORDERAPPPARTICLEADDON_H

/**
	@author Manu <manu@rustneversleeps>
*/

#include <vector>
#include <types.h>

class WithChannelOrderAppParticleAddon{
protected:
	std::vector<double> _channelOrderAPP;
public:
    WithChannelOrderAppParticleAddon(std::vector<double> channelOrderAPP);
    WithChannelOrderAppParticleAddon(int nChannelOrderAPP);

	WithChannelOrderAppParticleAddon(const WithChannelOrderAppParticleAddon& withChannelOrderAppParticleAddon);

	void SetChannelOrderAPP(double channelOrderAPP,int n) { _channelOrderAPP[n] = channelOrderAPP;}
	double GetChannelOrderAPP(int n) { return _channelOrderAPP[n];}

	int iMaxChannelOrderAPP();
};

#endif
