#include <iostream>
#include <math.h>
#include "Alphabet.h"
#include "exceptions.h"

// #define DEBUG

using namespace std;

Alphabet::Alphabet(int nBitsPorSimbolo,int longitudAlphabet,vector<vector<tBit> > secuenciasBits,vector<tSymbol> simbolos):_symbols(simbolos),_bitsSequences(secuenciasBits),_nBitsBySymbol(secuenciasBits[0].size()),_length(secuenciasBits.size())
{
    // si no coincide el numero de simbolos con el numero de secuencias de bits
    if(secuenciasBits.size()!=simbolos.size())
    {
			throw RuntimeException("Alphabet.cpp: el numero de secuencias de bits es distinto al de simbolos.");

    }

    computeMeanAndVariance();
}

Alphabet::Alphabet(vector<tSymbol> simbolos):_symbols(simbolos),_bitsSequences(simbolos.size(),vector<tBit>(0)),_nBitsBySymbol(0),_length(simbolos.size())
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

tSymbol Alphabet::operator [ ](vector<tBit> secuenciaBitsBuscada) const
{
    vector<vector<tBit> >::const_iterator iterator;
    iterator = find(_bitsSequences.begin(),_bitsSequences.end(),secuenciaBitsBuscada);
    if(iterator==_bitsSequences.end())
    {
			throw RuntimeException("Alphabet::operator[]: Esta secuencia de bits no forma parte del alfabeto.");
    }
	return _symbols[iterator - _bitsSequences.begin()];
}

vector<tBit> Alphabet::operator [ ](tSymbol simbolo) const
{
	vector<tSymbol>::const_iterator iterator;
	iterator = find(_symbols.begin(),_symbols.end(),simbolo);
	if(iterator==_symbols.end())
	{
		throw RuntimeException("Alphabet::operator[]: Este simbolo no forma parte del alfabeto.");
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
			throw RuntimeException("Alphabet::symbolsArray2int: Symbol not found.");
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
