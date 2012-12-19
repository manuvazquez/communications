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


#ifndef BASICTREENODE_H
#define BASICTREENODE_H

#include <vector>
#include <types.h>

class BasicTreeNode
{
protected:
	double _cost;
	uint _height,_id;
	std::vector<uint> _children;
	
public:	
    BasicTreeNode();
	BasicTreeNode(double cost);
	
	// defined so that the class is polymorphic
	virtual ~BasicTreeNode() {}
	
	double getCost() const { return _cost;}
	uint getHeight() const { return _height;}
	uint getID() const { return _id;}
	
	void setHeight(uint height) { _height = height;}
	void setID(uint id) { _id = id;}
	
	void addChild(uint id) { _children.push_back(id);}
	
	uint nChildren() const { return _children.size();}
};

#endif // BASICTREENODE_H
