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
        
        void computeMeanAndVariance();
    public:
        Alphabet(int nBitsPorSimbolo,int longitudAlphabet,std::vector<std::vector<tBit> > secuenciasBits,std::vector<tSymbol> simbolos);
        Alphabet(std::vector<tSymbol> simbolos);        
        int nBitsPerSymbol() const { return _nBitsBySymbol;}
        double variance() const { return _variance;}
        tSymbol operator [](std::vector<tBit> secuenciaBitsBuscada) const;
        tSymbol operator [](int index) const { return _symbols[index];}
        std::vector<tBit> operator [](tSymbol simbolo) const;
        int length() const { return _length;}
        void int2symbolsArray(int numero, std::vector<tSymbol> &res) const;
        int symbolsArray2int(std::vector<tSymbol> symbolsVector) const;
        tSymbol hardDecision(double softEstimation) const;
        tSymbol opposite(const tSymbol symbol) const { return -1.0*symbol;}
		bool doesItBelong(const tSymbol symbol) const;
		VectorXd int2eigenVector(int number, uint length) const;
		MatrixXd int2eigenMatrix(int number, uint rows, uint cols) const;
		
		/*!
		  it builds a new alphabet from this one with one additional symbol \ref symbol
		  \param symbol the symbol to be added
		  \return a new \ref Alphabet object
		*/
		Alphabet buildNewAlphabetByAddingSymbol(tSymbol symbol) const;
};
#endif

