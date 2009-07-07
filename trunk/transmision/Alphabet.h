#ifndef ALPHABET_H
#define ALPHABET_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <types.h>

class Alphabet
{
    private:
        std::vector<tSymbol> _symbols;
        std::vector<std::vector<tBit> > _bitsSequences;
        int _nBitsBySymbol,_length;
        double _mean,_variance;
    public:
        Alphabet(int nBitsPorSimbolo,int longitudAlphabet,std::vector<std::vector<tBit> > secuenciasBits,std::vector<tSymbol> simbolos);
        Alphabet(std::vector<tSymbol> simbolos);        
        int nBitsPerSymbol() const { return _nBitsBySymbol;}
        double variance() { return _variance;}
        tSymbol operator [](std::vector<tBit> secuenciaBitsBuscada);
        tSymbol operator [](int index) { return _symbols[index];}
        std::vector<tBit> operator [](tSymbol simbolo);
        int length() const { return _length;}
        void int2symbolsArray(int numero, std::vector<tSymbol> &res) const;
        int symbolsArray2int(std::vector<tSymbol> symbolsVector);
        tSymbol hardDecision(double softEstimation);
};
#endif

