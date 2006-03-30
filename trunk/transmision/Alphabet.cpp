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

    this->longitud = secuenciasBits.size();
    this->nBitsPorSimbolo = secuenciasBits[0].size();
    this->secuenciasBits = secuenciasBits;
    this->simbolos = simbolos;

    //se calcula la media y la varianza
    double media = 0;
    double mediaSimbolosCuadrado = 0;
    vector<tSymbol>::iterator iterador;
    for(iterador=simbolos.begin();iterador !=simbolos.end();iterador++)
    {
        media += (double) *iterador;
        mediaSimbolosCuadrado += ((double) *iterador)*((double)*iterador);
    }
    media /= longitud;
    mediaSimbolosCuadrado /= longitud;
    varianza = mediaSimbolosCuadrado - (media*media);
//     cout << "La media es " << media << " y la varianza " << varianza << endl;
}

tSymbol Alphabet::operator [ ](vector<tBit> secuenciaBitsBuscada)
{
    vector<vector<tBit> >::iterator iterador;
    iterador = find(secuenciasBits.begin(),secuenciasBits.end(),secuenciaBitsBuscada);
    if(iterador==secuenciasBits.end())
    {
			throw RuntimeException("Alphabet::operator[]: Esta secuencia de bits no forma parte del alfabeto.");
    }
//     cout << "Esta en la posicion " << iterador - secuenciasBits.begin() << endl;
	return simbolos[iterador - secuenciasBits.begin()];
}

vector<tBit> Alphabet::operator [ ](tSymbol simbolo)
{
	vector<tSymbol>::iterator iterador;
	iterador = find(simbolos.begin(),simbolos.end(),simbolo);
	if(iterador==simbolos.end())
	{
		throw RuntimeException("Alphabet::operator[]: Este simbolo no forma parte del alfabeto.");
	}
// 	cout << "Esta en la posicion" << (iterador - simbolos.begin()) << endl;
// 	printf("Secuencia Bits 0=%d,1=%d\n",secuenciasBits[iterador - simbolos.begin()][0],secuenciasBits[iterador - simbolos.begin()][1]);
	return secuenciasBits[iterador - simbolos.begin()];
}

void Alphabet::IntToArraySimbolos(int numero, vector<tSymbol> *res)
{
	int tamVector = res->size();

	if(numero>=pow((double)longitud,(double)tamVector))
		throw RuntimeException("El tamaño del vector es demasiado pequeño.");

// 	// se reserva espacio para un vector de bits tamaño "tamVector"
// 	vector<tBit> secuenciasBits(tamVector);
	
	int resto,i;

	i=1;
	do
	{
		resto = numero % longitud;
		cout << resto << "y" <<  simbolos[resto];
		(*res)[tamVector-i] =  simbolos[resto];
		numero /= longitud;
		i++;
	}while(numero!=0);

	for(;tamVector>=i;i++)
		(*res)[tamVector-i] = simbolos[0];
	cout << simbolos[0];
	
	cout << endl << "dentro" << endl;
	for(int j=0;j<(*res).size();j++)
		cout << (*res)[j];
	cout << endl;
}

