#include <iostream>
#include <cstdlib>
#include <vector>
#include "types.h"
#include <Alphabet.h>
#include <Bits.h>
#include <ARprocess.h>
#include <ARchannel.h>
#include <ChannelDependentNoise.h>
#include <Modulator.h>
#include <Demodulator.h>
#include <KalmanFilter.h>
#include <KalmanEstimator.h>
#include <ML_SMCAlgorithm.h>
#include <ResamplingCriterion.h>
#include <StatUtil.h>
#include <Util.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>
#include <mylapack.h>

using namespace std;
// using namespace la;
int main(int argc,char* argv[])
{
    int longitudAlphabet = 2;
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
//     tSymbol simbolos[] = {-3,-1,1,3};
    vector<vector<tBit> > secuenciasBits(4,vector<tBit>(2));
    secuenciasBits[0][0] = 0; secuenciasBits[0][1] = 0;
    secuenciasBits[1][0] = 0; secuenciasBits[1][1] = 1;
    secuenciasBits[2][0] = 1; secuenciasBits[2][1] = 0;
    secuenciasBits[3][0] = 1; secuenciasBits[3][1] = 1;
    vector<tSymbol> simbolos(4);
    simbolos[0] = -3; simbolos[1] = -1; simbolos[2] = 1; simbolos[3] = 3;
    Alphabet pam4(2,4,secuenciasBits,simbolos);

    vector<tBit> secuenciaAbuscar(2);
    secuenciaAbuscar[0] = 1; secuenciaAbuscar[1] = 1;
    tSymbol simboloDevuelto = pam4[secuenciaAbuscar];

	// -------------- PAM2 -------------------
//     vector<vector<tBit> > secuenciasBits(2,vector<tBit>(1));
	secuenciasBits = *(new vector<vector<tBit> >(2,vector<tBit>(1)));
    secuenciasBits[0][0] = 0; secuenciasBits[1][0] = 1;
//     vector<tSymbol> simbolos(2);
	simbolos = *(new vector<tSymbol>(2));
    simbolos[0] = -1; simbolos[1] = 1;
    Alphabet pam2(1,2,secuenciasBits,simbolos);


	vector<tBit> secuenciaDevuelta;
	secuenciaDevuelta = pam4[-1];
	cout << simboloDevuelto << endl;
// 	cout << "el tamaño de la secuencia es " << secuenciaDevuelta[0];
	secuenciaDevuelta.resize(6);
	vector<tSymbol> secuenciaSimbolos(8);
	pam4.IntToSymbolsArray(13,secuenciaSimbolos);
	cout << "Secuencia devuelta" << endl;
	for(int i=0;i<secuenciaSimbolos.size();i++)
		cout << secuenciaSimbolos[i];
	cout << endl;
	Bits bits(2,4);
	bits.Print();
	cout << "-----------" << endl;
	Bits bits2;

	tMatrix simbs = Modulator::Modulate(bits,pam2);

	cout << "Simbolos" << endl << simbs << "-------" << endl;

	Bits bitsDemodulados = Demodulator::Demodulate(simbs,pam2);

	cout << "Bits demodulados..." << endl;

	bitsDemodulados.Print();


	// ------------------------ Estimador de Kalman ------------------------------------
	int L=3,N=2,m=2,K=3;
	double channelMean=0.0,channelVariance=1.0,ARvariance=0.0001;
	vector<double> ARcoefficients(1);
	ARcoefficients[0] = 0.99999;

	Bits bitsTransmitir(N,K);

	tMatrix simbolosTransmitir = Modulator::Modulate(bitsTransmitir,pam2);

	tMatrix preambulo(N,m-1);
	preambulo = -1.0;
	simbolosTransmitir = Util::Append(preambulo,simbolosTransmitir);

	ARchannel canal(N,L,m,simbolosTransmitir.cols(),channelMean,channelVariance,ARcoefficients,ARvariance);
	cout << "longitud canal " << canal.Length() << " simbolos.co " << simbolosTransmitir.cols() << endl;
	for(int i=canal.Memory()-1;i<canal.Length();i++)
		cout << canal[i] << endl << "****************" << endl;

	ChannelDependentNoise ruido(canal);
	ruido.SetSNR(12,1);

	cout << "El ruido" << endl;
	ruido.Print();

	tMatrix observaciones = canal.Transmit(simbolosTransmitir,ruido);

	cout << "Las observaciones" << endl << observaciones;

	tMatrix mediaInicial(L,N*m);
	mediaInicial = 0.0;
	KalmanEstimator estimador(ARcoefficients[0],ARvariance,mediaInicial);

	// ----------------------------- Depuracion EstimadorKalman ----------------------------
// 	tMatrix matrizInicial(3,4);
// 	matrizInicial = 0.0;
// 	matrizInicial(0,0) = 0.1;matrizInicial(0,1) = 0.1342;matrizInicial(0,2) = 0.525;
// 	matrizInicial(1,0) = 1.1;matrizInicial(1,1) = -0.1342;matrizInicial(1,2) = 2.525;
// 	double coefAR = 0.99999;
// 	double varAR = 0.0001;
// 
// 
// 	KalmanEstimator estPrueba(coefAR,varAR,matrizInicial);
// 
// 	cout << "Los parametros" << endl << coefAR << "," << varAR << "," << endl << matrizInicial << endl;
// 
// 	tVector obs(3);
// 	obs(0) = 1.1; obs(1) = 1.3; obs(2) = 1.22;
// 	tMatrix simbolitos(2,2);
// 	simbolitos = 0.0;
// 	double vari = 0.52;
// 	tMatrix est = estPrueba.NextMatrix(obs,simbolitos,vari);
// 	cout << "los parametros de next: obs, simbolitos, vari" << endl << obs << endl << simbolitos << endl << vari << endl;
// 	cout << endl << "****mec mec******" << endl << est << endl << "**********" << endl;
// 	obs(0) = 3.1; obs(1) = 22.3; obs(2) = 0.22;
// 	simbolitos(0,1) = -1; simbolitos(0,0) = 13; simbolitos(1,0) = 4.0;
// 	est = estPrueba.NextMatrix(obs,simbolitos,1.0);
// 	cout << "los parametros de next: obs, simbolitos, vari" << endl << obs << endl << simbolitos << endl << vari << endl;
// 	cout << endl << "****mec mec******" << endl << est << endl << "**********" << endl;
// 	cout << "La verosimilitud es " << estPrueba.Likelihood(obs,simbolitos,1.0) << endl;
// 	simbolitos = 13.0;
// 	cout << "La verosimilitud es " << estPrueba.Likelihood(obs,simbolitos,1.0) << endl;
	// ------------------------------------------------------------------------------------

// 	KalmanEstimator estimador2 = *estimador.Clone();
// 
// 	tRange todasFilasSimbolos(0,N-1);
// 	for(int i=m-1;i<observaciones.cols();i++)
// 	{
// 		tRange rangoColumnas(i-m+1,i);
// 		tMatrix subMatrizSimbolos = simbolosTransmitir(todasFilasSimbolos,rangoColumnas);
// 		tMatrix est = estimador.NextMatrix(observaciones.col(i),subMatrizSimbolos,ruido.VarianceAt(i));
// 		cout << "Estimacion de Kalman (varianza es " << ruido.VarianceAt(i) << ")" << endl << est << endl;
// // 		cout << "El ruido" << endl << ruido
// 		cout << "Canal de verdad" << endl << canal[i] << endl << "-------------" << endl;
// 
// 		double verosimil = estimador.Likelihood(observaciones.col(i),subMatrizSimbolos,ruido.VarianceAt(i));
// 		cout << "La verosimilitud=" << verosimil << endl;
// 
// 		tMatrix simbolosChungos(N,m);
// 		simbolosChungos = -11;
// 		double verosimilChunga = estimador.Likelihood(observaciones.col(i),simbolosChungos,ruido.VarianceAt(i));
// 		cout << "La verosimilitud chunga=" << verosimilChunga << endl;		
// 	}

// ML_SMCAlgorithm::ML_SMCAlgorithm(string name, Alphabet alphabet, ChannelMatrixEstimator& channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingCriterion resamplingCriterion): SMCAlgorithm(name, alphabet, channelEstimator, preamble, smoothingLag, nParticles, resamplingCriterion)

// 	tMatrix preambulo(N,m-1);
// 	preambulo = -1.0;
	ML_SMCAlgorithm algoritmo("Detector suavizado optimo",pam2,estimador,preambulo,m-1,10,*(new ResamplingCriterion(0.9)));

// 	algoritmo.Run(observaciones,ruido.Variances());
	// --------------------------------------------------------------------------------------

	// ------------------------- Filtro de Kalman ------------------------------------------
// 	int longVectorAestimar = 4;
// 	int longVectorObservaciones = 3;
// 	double varEcEstado = 0.0001,varRuido = 0.2;
// 	double mediaVector=0.0,varVector=1.0;
// 	Random generador;
// // 	double* arrayNormal = generador.randnArray(longVectorAestimar*longVectorAestimar);
// // 	tMatrix R(arrayNormal,longVectorAestimar,longVectorAestimar);
// 	tMatrix R = LaGenMatDouble::eye(longVectorAestimar);
// 	R*=0.9999;
// 	double* arrayNormalF = generador.randnArray(longVectorObservaciones*longVectorAestimar);
// 	tMatrix F(arrayNormalF,longVectorObservaciones,longVectorAestimar);
// 	tMatrix covarianzaEcObservaciones = LaGenMatDouble::eye(longVectorObservaciones);
// 	covarianzaEcObservaciones*=varRuido;
// 	tMatrix covarianzaEcEstado = LaGenMatDouble::eye(longVectorAestimar);
// 	covarianzaEcEstado*=varEcEstado;
// 	tVector estado(generador.randnArray(longVectorAestimar,mediaVector,varVector),longVectorAestimar);
//
// 	tVector estadoInicialKalman(generador.randnArray(longVectorAestimar,mediaVector,varVector),longVectorAestimar);
// 	KalmanFilter kf(R,covarianzaEcEstado,estadoInicialKalman,LaGenMatDouble::eye(longVectorAestimar),longVectorObservaciones);
// 	tVector observacion(longVectorObservaciones);
// 	tVector nuevoEstado(longVectorAestimar);
// 	cout << "Estado al principio" << endl << estado << "---------" << endl;
// 	for(int i=0;i<30;i++)
// 	{
// 		Blas_Mat_Vec_Mult(R,estado,nuevoEstado);
// 		cout << "nuevoEstado" << endl << nuevoEstado << endl;
// 		for(int j=0;j<longVectorAestimar;j++)
// 			nuevoEstado(j) = nuevoEstado(j) + generador.randn()*varEcEstado;
// 		Blas_Mat_Vec_Mult(F,nuevoEstado,observacion);
// 		for(int j=0;j<longVectorObservaciones;j++)
// 			observacion(j) = observacion(j) + generador.randn()*varRuido;
// 		estado = nuevoEstado;
// 		kf.Step(F,observacion,covarianzaEcObservaciones);
// 		cout << "F,observacion,covarianzaEcObservaciones" << endl << F << endl << observacion << endl << covarianzaEcObservaciones << endl;
// 		cout << "El filtro de Kalman obtiene:" << endl << kf.PredictiveMean();
// 		cout << endl << "---------------" << endl;
// 		cout << estado;
// 	}
	//-------------------------------------------------------------------------------------------------------------


	// ------------------- De vector a Matrix ------------------------
// 	cout << "F es" << endl << F;
// 	tVector veci = Util::ToVector(F,rowwise);
// 	cout << "En forma de vector" << endl << veci << endl;
// 	cout << "De vuelta a la matrix" << endl << Util::ToMatrix(veci,rowwise,F.rows(),F.cols()) << endl;
	// -------------------------------------------------------------------

	// --------------------- Constructor copia de KalmanFilter --------------------------
// 	cout << kf.PredictiveCovariance() << endl;
// 	KalmanFilter kf2 = kf;
// 	kf2.Step(F,observacion,covarianzaEcObservaciones);
// 	cout << "------" << endl << kf.PredictiveCovariance() << endl;
// 	cout << "El modificado" << endl << kf2.PredictiveCovariance() << endl;
	// ---------------------------------------------------------------------------------






// 	// ----------- operaciones entre objetos Bits ------------------------
// 	Bits otrosBits(4,50);
// 	cout << "igual a los demodulados: "<< (bitsDemodulados==bits) << "a los otros: " << (bits==otrosBits) << endl;
//
// 	cout << "Diferencia: " << (bitsDemodulados-bits) << "a los otros: " << (bits-otrosBits) << endl;
// 	// -------------------------------------------------------------------

// 	bits2 = bits;
// 	bits2 = bits;
// 	bits2.Print();

// 	// ------------ Modulacion y demodulacion diferencial ------------
// 	Bits diffEncodBits = bits.DifferentialEncoding();
// 	diffEncodBits.Print();
// 	cout << "-----------" << endl;
// 	Bits diffDecodBits = diffEncodBits.DifferentialDecoding();
// 	diffDecodBits.Print();
// 	// --------------------------------------------------------------
//
// 	// -------------- Operaciones con matrices ----------------------
// 	tMatrix A(4,4); A = 1;
// 	tMatrix A2(2,4); A2 = 4.1;
// 	tMatrix A3(2,4);
// 	tMatrix B(4,4);
// 	tMatrix subA(2,2); subA = 0.1213;
// 	B = 2;
// 	Blas_Mat_Mat_Mult(A,B,A);
// 	cout << "Las matrices A y B" << endl << A << endl << B;
// 	tMatrix C = A*B;
// 	cout << A << endl << B << endl << C << endl;
// // 	Util::Add(A,A2,A3);
// // 	cout << A << endl << A2 << endl << A3 << endl;
// 	cout << "Apenddado" << endl << Util::Append(A,B) << endl;
// 	A(*(new tRange(1,2)),*(new tRange(2,3))).inject(subA);
// 	cout << "La nueva A" << endl << A << endl;
// 	// ----------------------------------------------------------------



//
// // 	Random r;
// // 	cout << r.randn() << endl;
//

//
// 	cout << "Ruido" << endl;
// 	ruido.Print();
// 	cout << "-----------" <<endl;
// 	ruido.SetSNR(3,1);
// 	ruido.Print();
// 	vector<double> varianzas = ruido.Variances();
// 	for(int i=0;i<ruido.Length();i++)
// 		cout << varianzas[i] << endl;
//
// 	cout << "una columna del ruido" << endl << ruido[3] << endl;
//
//
//
// 	tMatrix sub = matrizAleatoria(*(new tRange(1,2)),*(new tRange(1,2)));
// 	cout << "------- (is submatrix view" << sub.is_submatrixview() << ")" << endl << sub << endl;
//
// 	tMatrix otra(3,3);
// 	otra = 1;
// 	tVector v = otra.col(2);
//
// 	cout << "El vector es " << endl << v;
//
// 	tVector res(4);
// 	Blas_Mat_Vec_Mult(matrizAleatoria,v,res);
//
// 	cout << "matriz por vector" << endl << res;


	// ----------------- INVERSA -----------------------
// 	Random generador;
// 	double* arrayNormal = generador.randnArray(16);
// 	tMatrix matrizAleatoria(arrayNormal,4,4);
// 	cout << "Matriz aleatoria" << endl << matrizAleatoria << endl;
// 
// 	tLongIntVector pivotes(matrizAleatoria.size(0));
// 	LUFactorizeIP(matrizAleatoria,pivotes);
// 	double determinante = 1.0;
// 	for(int hola=0;hola<4;hola++)
// 		determinante *= matrizAleatoria(hola,hola);
// 	cout << "El determinante es " << determinante << endl;
// 	LaLUInverseIP(matrizAleatoria,pivotes);
// 
// 	cout << "La inversa" << endl << matrizAleatoria;
	// ------------------------------------------------


// 	// ------------------------------ proceso AR ---------------------------------------
// 	vector<double> coeficientes(1);
// 	coeficientes[0] = 0.99999;
// 	tMatrix matrizInicial(generador.randnArray(2*3),2,3);
// 	ARprocess procesoAR(matrizInicial,coeficientes,0.001);
// 	for(int i=0;i<200;i++)
// 	{
// 		cout << procesoAR.NextMatrix() << endl << "----------" << endl;
// 	}
// 	// ---------------------------------------------------------------------------------

	// --------------------------- Discrete_rnd --------------------------------------
// 	tVector probabilidades(3);
// 	probabilidades(0) = 0.6;
// 	probabilidades(1) = 1.2;
// 	probabilidades(2) = 0.6;
// 
// 	vector<int> indices = StatUtil::Discrete_rnd(10,probabilidades);
// 
// 	for(int i=0;i<indices.size();i++)
// 		cout << indices[i] << "---";
	// -------------------------------------------------------------------------------

    return 0;
}
