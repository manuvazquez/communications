#include <iostream>
#include <math.h>
#include "Alphabet.h"
#include "exceptions.h"

// #define DEBUG

using namespace std;

Alphabet::Alphabet(int nBitsPorSimbolo,int longitudAlphabet,vector<vector<tBit> > secuenciasBits,vector<tSymbol> simbolos)
{
    // si no coincide el numero de simbolos con el numero de secuencias de bits
    if(secuenciasBits.size()!=simbolos.size())
    {
			throw RuntimeException("Alphabet.cpp: el numero de secuencias de bits es distinto al de simbolos.");

    }

    this->_length = secuenciasBits.size();
    this->_nBitsBySymbol = secuenciasBits[0].size();
    this->_bitsSequences = secuenciasBits;
    this->_symbols = simbolos;

    //se calcula la media y la varianza
    double _mean = 0;
    double mediaSimbolosCuadrado = 0;
    vector<tSymbol>::iterator iterador;
    for(iterador=simbolos.begin();iterador !=simbolos.end();iterador++)
    {
        _mean += (double) *iterador;
        mediaSimbolosCuadrado += ((double) *iterador)*((double)*iterador);
    }
    _mean /= _length;
    mediaSimbolosCuadrado /= _length;
    _variance = mediaSimbolosCuadrado - (_mean*_mean);
//     cout << "Mean is " << _mean << " and variance " << _variance << endl;
}

tSymbol Alphabet::operator [ ](vector<tBit> secuenciaBitsBuscada)
{
    vector<vector<tBit> >::iterator iterador;
    iterador = find(_bitsSequences.begin(),_bitsSequences.end(),secuenciaBitsBuscada);
    if(iterador==_bitsSequences.end())
    {
			throw RuntimeException("Alphabet::operator[]: Esta secuencia de bits no forma parte del alfabeto.");
    }
//     cout << "Esta en la posicion " << iterador - _bitsSequences.begin() << endl;
	return _symbols[iterador - _bitsSequences.begin()];
}

vector<tBit> Alphabet::operator [ ](tSymbol simbolo)
{
	vector<tSymbol>::iterator iterador;
	iterador = find(_symbols.begin(),_symbols.end(),simbolo);
	if(iterador==_symbols.end())
	{
		throw RuntimeException("Alphabet::operator[]: Este simbolo no forma parte del alfabeto.");
	}
	return _bitsSequences[iterador - _symbols.begin()];
}

void Alphabet::IntToSymbolsArray(int numero, vector<tSymbol> &res)
{
	int tamVector = res.size();

	#ifdef DEBUG
		cout << "tamVector: " << tamVector << endl;
	#endif

	if(numero>=pow((double)_length,(double)tamVector))
		throw RuntimeException("Alphabet::IntToSymbolsArray: vector size is smaller than needed.");

// 	// se reserva espacio para un vector de bits tama�o "tamVector"
// 	vector<tBit> _bitsSequences(tamVector);

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

int Alphabet::SymbolsArrayToInt(vector<tSymbol> symbolsVector)
{
	int size = symbolsVector.size();

	#ifdef DEBUG
		cout << "symbolsVector.size(): " << symbolsVector.size() << endl;
	#endif

	int res = 0, base = 1;
	vector<tSymbol>::iterator iterator;
	for(int i=size-1;i>=0;i--)
	{
		#ifdef DEBUG
			cout << "Antes del for" << endl;
		#endif
		iterator = find(_symbols.begin(),_symbols.end(),symbolsVector.at(i));
		if(iterator==_symbols.end())
		{
			throw RuntimeException("Alphabet::SymbolsArrayToInt: Symbol not found.");
		}
		res += base*(iterator - _symbols.begin());
		base *= _length;
	}
	return res;
}

tSymbol Alphabet::HardDecision(double softEstimation)
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
