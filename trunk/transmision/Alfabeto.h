#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <vector>
#include "types.h"

using namespace std;
class Alfabeto
{
	private:
		vector<tSymbol> simbolos;
		vector<vector<tBit> > secuenciasBits;
		int nBitsPorSimbolo,longitud;
		double media,varianza;
	public:
		Alfabeto(int nBitsPorSimbolo,int longitudAlfabeto,vector<vector<tBit> > secuenciasBits,vector<tSymbol> simbolos);
		int NbitsPorSimbolo() { return nBitsPorSimbolo;}
		double Varianza() { return varianza;}
		tSymbol operator [](vector<tBit> secuenciaBitsBuscada);
		vector<tBit> operator [](tSymbol simbolo);
		int Longitud() { return longitud;}
		void IntToArraySimbolos(int numero, vector<tSymbol> *res);
};
#endif

