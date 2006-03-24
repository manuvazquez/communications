#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <vector>
#include "tipos.h"

using namespace std;
class Alfabeto
{
	private:
		vector<tSimbolo> simbolos;
		vector<vector<tBit> > secuenciasBits;
		int nBitsPorSimbolo,longitud;
		double media,varianza;
	public:
		Alfabeto(int nBitsPorSimbolo,int longitudAlfabeto,vector<vector<tBit> > secuenciasBits,vector<tSimbolo> simbolos);
		int NbitsPorSimbolo();
		double Varianza();
		tSimbolo operator [](vector<tBit> secuenciaBitsBuscada);
		vector<tBit> operator [](tSimbolo simbolo);
		int Longitud();
		void IntToArraySimbolos(int numero, vector<tSimbolo> *res);
};
#endif

