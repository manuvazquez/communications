#include <iostream>
#include <cstdlib>
#include <vector>
#include "tipos.h"
#include "Alfabeto.h"
#include <Bits.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas3pp.h>

using namespace std;
// using namespace la;
int main(int argc,char* argv[])
{
    int longitudAlfabeto = 2;
//     LaGenMatDouble matriz(4,2);
//     LaGenMatDouble matriz2(2,3);
//     LaGenMatDouble res(4,3);
//     matriz = 3;
//     matriz2 = 1;
//     matriz2(1,1) = 3;
//     Blas_Mat_Mat_Mult(matriz,matriz2,res);
// //     matriz2(2,2) = 18;
//     cout << "Hola, holita\n" << matriz << matriz2 << "a*b\n" << res;

//     cout << "El primer argumento " << argv[1] << "\n";
//     tBit secuenciasBits[][2] = {{0,0},{0,1},{1,0},{1,1}};
//     tSimbolo simbolos[] = {-3,-1,1,3};
    vector<vector<tBit> > secuenciasBits(4,vector<tBit>(2));
    secuenciasBits[0][0] = 0; secuenciasBits[0][1] = 0;
    secuenciasBits[1][0] = 0; secuenciasBits[1][1] = 1;
    secuenciasBits[2][0] = 1; secuenciasBits[2][1] = 0;
    secuenciasBits[3][0] = 1; secuenciasBits[3][1] = 1;
    vector<tSimbolo> simbolos(4);
    simbolos[0] = -3; simbolos[1] = -1; simbolos[2] = 1; simbolos[3] = 3;
    Alfabeto pam4(2,4,secuenciasBits,simbolos);
    vector<tBit> secuenciaAbuscar(2);
    secuenciaAbuscar[0] = 1; secuenciaAbuscar[1] = 1;
    tSimbolo simboloDevuelto = pam4[secuenciaAbuscar];

	vector<tBit> secuenciaDevuelta;
	secuenciaDevuelta = pam4[-1];
	cout << simboloDevuelto << endl;
// 	cout << "el tamaño de la secuencia es " << secuenciaDevuelta[0];
	secuenciaDevuelta.resize(6);
	vector<tSimbolo> secuenciaSimbolos(8);
	pam4.IntToArraySimbolos(13,&secuenciaSimbolos);
	cout << "Secuencia devuelta" << endl;
	for(int i=0;i<secuenciaSimbolos.size();i++)
		cout << secuenciaSimbolos[i];
	cout << endl;
	Bits bits(4,2);
	bits.Print();
	cout << "-----------" << endl;
	Bits bits2;

// 	bits2 = bits;
// 	bits2 = bits;
// 	bits2.Print();

	Bits diffEncodBits = bits.DifferentialEncoding();
	diffEncodBits.Print();
	cout << "-----------" << endl;
	Bits diffDecodBits = diffEncodBits.DifferentialDecoding();
	diffDecodBits.Print();
// 	cin >> caracter;
    return 0;
}
