#ifndef ALPHABET_H
#define ALPHABET_H

#include <iostream>
#include <vector>
#include <types.h>

using namespace std;
class Alphabet
{
	private:
		vector<tSymbol> _symbols;
		vector<vector<tBit> > _bitsSequences;
		int _nBitsBySymbol,_length;
		double _mean,_variance;
	public:
		Alphabet(int nBitsPorSimbolo,int longitudAlphabet,vector<vector<tBit> > secuenciasBits,vector<tSymbol> simbolos);
		int NbitsBySymbol() { return _nBitsBySymbol;}
		double Variance() { return _variance;}
		tSymbol operator [](vector<tBit> secuenciaBitsBuscada);
		tSymbol operator [](int index) { return _symbols[index];}
		vector<tBit> operator [](tSymbol simbolo);
		int Length() { return _length;}
		void IntToSymbolsArray(int numero, vector<tSymbol> &res);
};
#endif

