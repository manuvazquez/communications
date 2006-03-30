#include <iostream>
#include <math.h>
#include "Alphabet.h"
#include "exceptions.h"

using namespace std;

Alphabet::Alphabet(int nBitsPorSimbolo,int longitudAlphabet,vector<vector<tBit> > secuenciasBits,vector<tSymbol> simbolos)
{
    cout << "La longitud del alfabeto es: " << secuenciasBits.size() << endl;
    cout << "y el numero de bits por simbolo " << secuenciasBits[0].size() << endl;

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
// 	cout << "Esta en la posicion" << (iterador - _symbols.begin()) << endl;
// 	printf("Secuencia Bits 0=%d,1=%d\n",_bitsSequences[iterador - _symbols.begin()][0],_bitsSequences[iterador - _symbols.begin()][1]);
	return _bitsSequences[iterador - _symbols.begin()];
}

void Alphabet::IntToSymbolsArray(int numero, vector<tSymbol> *res)
{
	int tamVector = res->size();

	if(numero>=pow((double)_length,(double)tamVector))
		throw RuntimeException("El tamaño del vector es demasiado pequeño.");

// 	// se reserva espacio para un vector de bits tamaño "tamVector"
// 	vector<tBit> _bitsSequences(tamVector);
	
	int resto,i;

	i=1;
	do
	{
		resto = numero % _length;
		cout << resto << "y" <<  _symbols[resto];
		(*res)[tamVector-i] =  _symbols[resto];
		numero /= _length;
		i++;
	}while(numero!=0);

	for(;tamVector>=i;i++)
		(*res)[tamVector-i] = _symbols[0];
	cout << _symbols[0];
	
	cout << endl << "dentro" << endl;
	for(int j=0;j<(*res).size();j++)
		cout << (*res)[j];
	cout << endl;
}

