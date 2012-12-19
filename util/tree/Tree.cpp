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


#include "Tree.h"

#include <iostream>
#include <typeinfo>

Tree::~Tree()
{
	for(uint i=0;i<_nodes.size();i++)
		delete _nodes[i];
}

void Tree::setRoot(BasicTreeNode *node)
{
	_nodes.push_back(node);
}

void Tree::addNode(BasicTreeNode* child, BasicTreeNode* parent)
{
	uint id = _nodes.size();
	
	child->setHeight(parent->getHeight()+1);
	child->setID(id);
	_nodes.push_back(child);
	
	parent->addChild(id);
}

BasicTreeNode* Tree::bestLeaf()
{
	assert(_nodes.size()>0);
    uint iBest;
    double bestCost = 0.0;

	uint i=0;
	
	// the first leaf node is searched for
	while(i<_nodes.size() && _nodes[i]->nChildren()!=0)
		i++;
	
	iBest = i;
	bestCost = _nodes[i]->getCost();
	
    for(i=i+1;i<_nodes.size();i++)
    {
        // if it isn't a leaf node...
        if(_nodes[i]->nChildren()!=0)
            continue;

        if(_nodes[i]->getCost() < bestCost)
        {
            iBest = i;
            bestCost = _nodes[i]->getCost();
        }
    }
        
    return _nodes[iBest];
}