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
#include "WithChannelOrderAppParticleAddon.h"

WithChannelOrderAppParticleAddon::WithChannelOrderAppParticleAddon(std::vector<double> channelOrderAPP):_channelOrderAPP(channelOrderAPP)
{
}

WithChannelOrderAppParticleAddon::WithChannelOrderAppParticleAddon(int nChannelOrderAPP):_channelOrderAPP(nChannelOrderAPP)
{
	for(int i=0;i<_channelOrderAPP.size();i++)
		_channelOrderAPP[i] = 1.0/(double)_channelOrderAPP.size();
}

WithChannelOrderAppParticleAddon::WithChannelOrderAppParticleAddon(const WithChannelOrderAppParticleAddon& withChannelOrderAppParticleAddon):_channelOrderAPP(withChannelOrderAppParticleAddon._channelOrderAPP)
{
}
