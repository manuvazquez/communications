#include <iostream>
#include <math.h>
#include "Alphabet.h"
#include "exceptions.h"

// #define DEBUG

using namespace std;

Alphabet::Alphabet(vector<vector<tBit> > bitsSequences,vector<tSymbol> symbols):_symbols(symbols),_bitsSequences(bitsSequences),_nBitsPerSymbol(bitsSequences[0].size()),_length(bitsSequences.size())
{
    // if the number of bits sequences and that of the symbols don't match...
    if(bitsSequences.size()!=symbols.size())
	  throw RuntimeException("Alphabet::Alphabet: number of bit sequences is different from that of symbols.");

    computeMeanAndVariance();
}

Alphabet::Alphabet(vector<tSymbol> symbols):_symbols(symbols),_bitsSequences(symbols.size(),vector<tBit>(0)),_nBitsPerSymbol(0),_length(symbols.size())
{
    computeMeanAndVariance();
}

void Alphabet::computeMeanAndVariance()
{
    _mean = 0.0;
    double squaredSymbolsMean = 0.0;
    vector<tSymbol>::iterator iterator;
    for(iterator=_symbols.begin();iterator !=_symbols.end();iterator++)
    {
        _mean += (double) *iterator;
        squaredSymbolsMean += ((double) *iterator)*((double)*iterator);
    }
    _mean /= _length;
    squaredSymbolsMean /= _length;
    _variance = squaredSymbolsMean - (_mean*_mean);
}

tSymbol Alphabet::operator [ ](vector<tBit> searchedBitsSequence) const
{
    vector<vector<tBit> >::const_iterator iterator;
    iterator = find(_bitsSequences.begin(),_bitsSequences.end(),searchedBitsSequence);
    if(iterator==_bitsSequences.end())
    {
			throw RuntimeException("Alphabet::operator[]: this sequence of bits doesn't belong to the alphabet.");
    }
	return _symbols[iterator - _bitsSequences.begin()];
}

vector<tBit> Alphabet::operator [ ](tSymbol symbol) const
{
	vector<tSymbol>::const_iterator iterator;
	iterator = find(_symbols.begin(),_symbols.end(),symbol);
	if(iterator==_symbols.end())
	{
		throw RuntimeException("Alphabet::operator[]: this symbol doesn't belong to the alphabet.");
	}
	return _bitsSequences[iterator - _symbols.begin()];
}

void Alphabet::int2symbolsArray(int numero, vector<tSymbol> &res) const
{
	int tamVector = res.size();

	if(numero>=pow((double)_length,(double)tamVector))
		throw RuntimeException("Alphabet::int2symbolsArray: vector size is smaller than needed.");

	int resto,i;

	i=1;
	do
	{
		resto = numero % _length;
		res[tamVector-i] =  _symbols[resto];
		numero /= _length;
		i++;
	}while(numero!=0);

	for(;tamVector>=i;i++)
		res[tamVector-i] = _symbols[0];
}

int Alphabet::symbolsArray2int(vector<tSymbol> symbolsVector) const
{
	int size = symbolsVector.size();

	int res = 0, base = 1;
	vector<tSymbol>::const_iterator iterator;
	for(int i=size-1;i>=0;i--)
	{
		iterator = find(_symbols.begin(),_symbols.end(),symbolsVector.at(i));
		if(iterator==_symbols.end())
		{
			throw RuntimeException("Alphabet::symbolsArray2int: symbol not found.");
		}
		res += base*(iterator - _symbols.begin());
		base *= _length;
	}
	return res;
}

tSymbol Alphabet::hardDecision(double softEstimation) const
{
	double distance;

	double minDistance = fabs(softEstimation - _symbols[0]);
	int iMin = 0;

	for(int i=1;i<_length;i++)
	{
		distance = fabs(softEstimation - _symbols[i]);
		if(distance<minDistance)
		{
			minDistance = distance;
			iMin = i;
		}
	}
	return _symbols[iMin];
}

bool Alphabet::doesItBelong(const tSymbol symbol) const
{
	vector<tSymbol>::const_iterator iterator;
	
	iterator = find(_symbols.begin(),_symbols.end(),symbol);
	
	if(iterator==_symbols.end())
	  return false;
	else
	  return true;
}

VectorXd Alphabet::int2eigenVector(int number, uint length) const
{
	if(number>=pow((double)_length,(double)length))
		throw RuntimeException("Alphabet::int2eigenVector: vector size is smaller than needed.");

	VectorXd res = VectorXd::Constant(length,_symbols[0]);
	
	uint remainder,i;

	i=1;
	do
	{
		remainder = number % _length;
		res(length-i) =  _symbols[remainder];
		number /= _length;
		i++;
	}while(number!=0);

	return res;
}

MatrixXd Alphabet::int2eigenMatrix(int number, uint rows, uint cols) const
{
	uint length = rows*cols;
	
	if(number>=pow((double)_length,(double)(length)))
		throw RuntimeException("Alphabet::int2eigenVector: vector size is smaller than needed.");

	MatrixXd res = MatrixXd::Constant(rows,cols,_symbols[0]);
	
	uint remainder;

	length--;
	do
	{
		remainder = number % _length;
		res(length % rows,length / rows) = _symbols[remainder];
		number /= _length;
		length--;
	}while(number!=0);

	return res;
}

Alphabet Alphabet::buildNewAlphabetByAddingSymbol(tSymbol symbol) const
{
  vector<tSymbol> newSymbols = _symbols;
  
  newSymbols.push_back(symbol);
  
  return Alphabet(newSymbols);
}