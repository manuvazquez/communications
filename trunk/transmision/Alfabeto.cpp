#include <iostream>
#include <math.h>
#include "Alfabeto.h"
#include "excepcionesTransmision.h"

using namespace std;

Alfabeto::Alfabeto(int nBitsPorSimbolo,int longitudAlfabeto,vector<vector<tBit> > secuenciasBits,vector<tSimbolo> simbolos)
{
    cout << "La longitud del alfabeto es: " << secuenciasBits.size() << endl;
    cout << "y el numero de bits por simbolo " << secuenciasBits[0].size() << endl;

    // si no coincide el numero de simbolos con el numero de secuencias de bits
    if(secuenciasBits.size()!=simbolos.size())
    {
			throw RuntimeException("Alfabeto.cpp: el numero de secuencias de bits es distinto al de simbolos.");
			
    }

    this->longitud = secuenciasBits.size();
    this->nBitsPorSimbolo = secuenciasBits[0].size();
    this->secuenciasBits = secuenciasBits;
    this->simbolos = simbolos;

    //se calcula la media y la varianza
    double media = 0;
    double mediaSimbolosCuadrado = 0;
    vector<tSimbolo>::iterator iterador;
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

int Alfabeto::NbitsPorSimbolo()
{
    return nBitsPorSimbolo;
}

double Alfabeto::Varianza()
{
    return varianza;
}

tSimbolo Alfabeto::operator [ ](vector<tBit> secuenciaBitsBuscada)
{
    vector<vector<tBit> >::iterator iterador;
    iterador = find(secuenciasBits.begin(),secuenciasBits.end(),secuenciaBitsBuscada);
    if(iterador==secuenciasBits.end())
    {
			throw RuntimeException("Alfabeto::operator[]: Esta secuencia de bits no forma parte del alfabeto.");
    }
    cout << "Esta en la posicion " << iterador - secuenciasBits.begin() << endl;
	return simbolos[iterador - secuenciasBits.begin()];
}

vector<tBit> Alfabeto::operator [ ](tSimbolo simbolo)
{
	vector<tSimbolo>::iterator iterador;
	iterador = find(simbolos.begin(),simbolos.end(),simbolo);
	if(iterador==simbolos.end())
	{
		throw RuntimeException("Alfabeto::operator[]: Este simbolo no forma parte del alfabeto.");
	}
	cout << "Esta en la posicion" << (iterador - simbolos.begin()) << endl;
	printf("Secuencia Bits 0=%d,1=%d\n",secuenciasBits[iterador - simbolos.begin()][0],secuenciasBits[iterador - simbolos.begin()][1]);
	return secuenciasBits[iterador - simbolos.begin()];
}

int Alfabeto::Longitud()
{
	return longitud;
}

void Alfabeto::IntToArraySimbolos(int numero, vector<tSimbolo> *res)
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

// public sbyte[] IntToArraySimbolos(int numero, int tamVector)
// {
// 	if(numero>=Math.Pow(longAlfabeto,tamVector))
// 		throw new System.ApplicationException("El tamaÃ±o del vector es demasiado pequeÃ±o para representar ese entero.");
// 	else
// 	{
// 		sbyte[] array = new sbyte[tamVector];
// 		int resto,i;
//
// 		i = 1;
// 		do
// 		{
// 				resto = numero % longAlfabeto;
// 				array[tamVector-i]=simbolosCorrespondientes[resto];
// 				numero = numero / longAlfabeto;
// 				i++;
// 		}while(numero!=0);
// 		for(;tamVector>=i;i++)
// 			array[tamVector-i] = simbolosCorrespondientes[0];
// 		return array;
// 	}
// }


// public void IncrementarArraySimbolos(sbyte[] vectorSimbolos)
// {
// 	int j;
// 	int nSimbolos = vectorSimbolos.Length;
// 	sbyte simboloAbuscar;
//
// 	for(int i=nSimbolos-1;i>-1;i--)
// 	{
// 		simboloAbuscar = vectorSimbolos[i];
//
// 		// se busca el simbolo
// 		j = 0;
// 		while(simbolosCorrespondientes[j]!=simboloAbuscar)
// 			j++;
//
// 		// si no es el ultimo
// 		if(j<(longAlfabeto-1))
// 		{
// 			vectorSimbolos[i] = simbolosCorrespondientes[j+1];
// 			return;
// 		}
// 		else
// 			vectorSimbolos[i] = simbolosCorrespondientes[0];
// 	}
// }
//
// public sbyte[] PrimerArraySimbolos(int longitud)
// {
// 	sbyte[] arraySimbolos = new sbyte[longitud];
// 	for(int i=0;i<longitud;i++)
// 		arraySimbolos[i] = simbolosCorrespondientes[0];
// 	return arraySimbolos;
// }
//
