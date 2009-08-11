#ifndef TIPOS_H
#define TIPOS_H

#include <lapackpp/gmd.h>
#include <lapackpp/lavd.h>
#include <lapackpp/laindex.h>
#include <lapackpp/lavli.h>

#include <Eigen/Core>

typedef double tSymbol;
typedef unsigned short int tBit;
// typedef char tBit;
typedef LaGenMatDouble tMatrix;
typedef LaVectorDouble tVector;
typedef LaVectorLongInt tLongIntVector;
typedef LaIndex tRange;

typedef unsigned int uint;

USING_PART_OF_NAMESPACE_EIGEN

#endif
