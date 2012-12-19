/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2012  Manu <manuavazquez@gmail.com>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#ifndef KNOWNFLATCHANNELOPTIMALALGORITHMNODE_H
#define KNOWNFLATCHANNELOPTIMALALGORITHMNODE_H

#include <BasicTreeNode.h>


class KnownFlatChannelOptimalAlgorithmNode : public BasicTreeNode
{
protected:
	VectorXd _symbolsVector;
public:
    KnownFlatChannelOptimalAlgorithmNode(uint nSymbols);
    KnownFlatChannelOptimalAlgorithmNode(double cost, VectorXd symbolsVector);
	
	VectorXd getSymbolsVector() const { return _symbolsVector;}
};

#endif // KNOWNFLATCHANNELOPTIMALALGORITHMNODE_H
