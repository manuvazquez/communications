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
#include "PSPPath.h"

// #define DEBUG

PSPPath::PSPPath(): ViterbiPath(),_estimatedChannelMatrices(NULL)
{
}


PSPPath::PSPPath(int nTimeInstants,double cost, tMatrix initialSequence, std::vector<std::vector<tMatrix> > initialChannelMatrices, std::vector<ChannelMatrixEstimator *> channelMatrixEstimators): ViterbiPath(nTimeInstants, cost, initialSequence), _channelMatrixEstimators(channelMatrixEstimators.size())
{
	if(initialChannelMatrices.size()!=channelMatrixEstimators.size())
		throw RuntimeException("PSPPath::PSPPath: channel order implied by the length of the \"initialChannelMatrices\" vector is not equal to that implied by \"channelMatrixEstimators\".");

	if(initialChannelMatrices[0].size()>initialSequence.cols())
		throw RuntimeException("PSPPath::PSPPath: number of received detected symbol vectors is less than number of received detected channel matrices.");

	_estimatedChannelMatrices = new tMatrix*[channelMatrixEstimators.size()];

	for(uint iChannelMatrixEstimator=0;iChannelMatrixEstimator<channelMatrixEstimators.size();iChannelMatrixEstimator++)
	{
		_channelMatrixEstimators[iChannelMatrixEstimator] = channelMatrixEstimators[iChannelMatrixEstimator]->Clone();

        _estimatedChannelMatrices[iChannelMatrixEstimator] = new tMatrix[_nTimeInstants];
		for(uint i=0;i<initialChannelMatrices[iChannelMatrixEstimator].size();i++)
			_estimatedChannelMatrices[iChannelMatrixEstimator][initialSequence.cols()-1-i] = initialChannelMatrices[iChannelMatrixEstimator][i];
	}
}

PSPPath::PSPPath(const PSPPath &path):ViterbiPath(path),_channelMatrixEstimators(path._channelMatrixEstimators.size()),_estimatedChannelMatrices(new tMatrix*[path._channelMatrixEstimators.size()])
{
	for(uint iChannelMatrixEstimator=0;iChannelMatrixEstimator<_channelMatrixEstimators.size();iChannelMatrixEstimator++)
	{
		_channelMatrixEstimators[iChannelMatrixEstimator] = path._channelMatrixEstimators[iChannelMatrixEstimator]->Clone();
        _estimatedChannelMatrices[iChannelMatrixEstimator] = new tMatrix[_nTimeInstants];
			for(int i=0;i<_nTimeInstants;i++)
				_estimatedChannelMatrices[iChannelMatrixEstimator][i] = path._estimatedChannelMatrices[iChannelMatrixEstimator][i];
	}
}


PSPPath::~PSPPath()
{
    for(uint i=0;i<_channelMatrixEstimators.size();i++)
    {
        delete _channelMatrixEstimators[i];
        delete[] _estimatedChannelMatrices[i];
    }
    delete[] _estimatedChannelMatrices;
}


void PSPPath::Clean()
{
    ViterbiPath::Clean();
	for(uint i=0;i<_channelMatrixEstimators.size();i++)
	{
		delete _channelMatrixEstimators[i];
		_channelMatrixEstimators[i] = NULL;

		delete[] _estimatedChannelMatrices[i];
		_estimatedChannelMatrices[i] = NULL;

	}
//     delete[] _estimatedChannelMatrices;
//     _estimatedChannelMatrices = NULL;
}

void PSPPath::Print() const
{
    ViterbiPath::Print();
    cout << "number of channel matrix estimators: " << _channelMatrixEstimators.size() << endl;
	for(uint iChannelMatrixEstimator=0;iChannelMatrixEstimator<_channelMatrixEstimators.size();iChannelMatrixEstimator++)
	{
		cout << "channel order index: " << iChannelMatrixEstimator << endl << _estimatedChannelMatrices[iChannelMatrixEstimator][_detectedSequence->cols()-1] << endl;
	}
}

void PSPPath::Update(const PSPPath& path, tVector newSymbolVector, double newCost, std::vector<ChannelMatrixEstimator *> newChannelMatrixEstimators/*,const std::vector<tMatrix> &newChannelMatrices*/)
{
	if(newChannelMatrixEstimators.size()!=path._channelMatrixEstimators.size())
		throw RuntimeException("PSPPath::Update: the number of ChannelMatrixEstimator's fo the source path object and the number of the received ones differ.");

    #ifdef DEBUG
    	cout << "Antes de llamar al Update de ViterbiPath" << endl;
    	cout << "newSymbolVector es " << endl << newSymbolVector;
    	cout << "newCost es " << newCost << endl;
    #endif

    ViterbiPath::Update(path, newSymbolVector, newCost);

    #ifdef DEBUG
    	cout << "LLamado al Update de ViterbiPath" << endl;
    #endif

	// if this object does not have the proper number of ChannelMatrixEstimator's
	if(_channelMatrixEstimators.size()==0)
	{
		_channelMatrixEstimators.resize(path._channelMatrixEstimators.size(),NULL);
	}
	if(_estimatedChannelMatrices==NULL)
	{
		_estimatedChannelMatrices = new tMatrix*[_channelMatrixEstimators.size()];
		for(uint iChannelMatrixEstimator=0;iChannelMatrixEstimator<_channelMatrixEstimators.size();iChannelMatrixEstimator++)
			_estimatedChannelMatrices[iChannelMatrixEstimator] = NULL;
	}

	for(uint iChannelMatrixEstimator=0;iChannelMatrixEstimator<path._channelMatrixEstimators.size();iChannelMatrixEstimator++)
	{
		delete _channelMatrixEstimators[iChannelMatrixEstimator];
		_channelMatrixEstimators[iChannelMatrixEstimator] = newChannelMatrixEstimators[iChannelMatrixEstimator];

		delete[] _estimatedChannelMatrices[iChannelMatrixEstimator];
		_estimatedChannelMatrices[iChannelMatrixEstimator] = new tMatrix[_nTimeInstants];

		for(int i=0;i<_nTimeInstants;i++)
			_estimatedChannelMatrices[iChannelMatrixEstimator][i] = path._estimatedChannelMatrices[iChannelMatrixEstimator][i];

		// the new matrix is added at the right index based on the last detected vector
		_estimatedChannelMatrices[iChannelMatrixEstimator][_detectedSequence->cols()-1] = newChannelMatrixEstimators[iChannelMatrixEstimator]->LastEstimatedChannelMatrix();
	}
}

void PSPPath::operator=(const PSPPath &path)
{
	ViterbiPath::operator =(path);

	// this Path has not been intialized
	if(_estimatedChannelMatrices==NULL)
	{
		_estimatedChannelMatrices = new tMatrix*[path._channelMatrixEstimators.size()];
		_channelMatrixEstimators.resize(path._channelMatrixEstimators.size());
		for(uint iChannelMatrixEstimator=0;iChannelMatrixEstimator<path._channelMatrixEstimators.size();iChannelMatrixEstimator++)
		{
			_channelMatrixEstimators[iChannelMatrixEstimator] = path._channelMatrixEstimators[iChannelMatrixEstimator]->Clone();

			_estimatedChannelMatrices[iChannelMatrixEstimator] = new tMatrix[_nTimeInstants];
			for(int i=0;i<_nTimeInstants;i++)
				_estimatedChannelMatrices[iChannelMatrixEstimator][i] = path._estimatedChannelMatrices[iChannelMatrixEstimator][i];
		}
	}else
	{
		if(_channelMatrixEstimators.size()!=path._channelMatrixEstimators.size())
			throw RuntimeException("PSPPath::operator=: Paths objects being equaled have different number of channel matrix estimators.");
		for(uint iChannelMatrixEstimator=0;iChannelMatrixEstimator<path._channelMatrixEstimators.size();iChannelMatrixEstimator++)
		{
			delete _channelMatrixEstimators[iChannelMatrixEstimator];
			_channelMatrixEstimators[iChannelMatrixEstimator] = path._channelMatrixEstimators[iChannelMatrixEstimator]->Clone();

			delete[] _estimatedChannelMatrices[iChannelMatrixEstimator];
			_estimatedChannelMatrices[iChannelMatrixEstimator] = new tMatrix[_nTimeInstants];
			for(int i=0;i<_nTimeInstants;i++)
				_estimatedChannelMatrices[iChannelMatrixEstimator][i] = path._estimatedChannelMatrices[iChannelMatrixEstimator][i];
		}
    }
}
