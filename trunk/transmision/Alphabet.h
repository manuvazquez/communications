#ifndef ALPHABET_H
#define ALPHABET_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <types.h>

// using namespace std;
class Alphabet
{
	private:
		std::vector<tSymbol> _symbols;
		std::vector<std::vector<tBit> > _bitsSequences;
		int _nBitsBySymbol,_length;
		double _mean,_variance;
	public:
		Alphabet(int nBitsPorSimbolo,int longitudAlphabet,std::vector<std::vector<tBit> > secuenciasBits,std::vector<tSymbol> simbolos);
		int NbitsBySymbol() { return _nBitsBySymbol;}
		double Variance() { return _variance;}
		tSymbol operator [](std::vector<tBit> secuenciaBitsBuscada);
		tSymbol operator [](int index) { return _symbols[index];}
		std::vector<tBit> operator [](tSymbol simbolo);
		int Length() const { return _length;}
		void IntToSymbolsArray(int numero, std::vector<tSymbol> &res) const;
		int SymbolsArrayToInt(std::vector<tSymbol> symbolsVector);
		tSymbol HardDecision(double softEstimation);
};
#endif

