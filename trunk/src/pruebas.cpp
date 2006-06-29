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
#include <RLSEstimator.h>
#include <LMSEstimator.h>
#include <RMMSEDetector.h>
#include <ML_SMCAlgorithm.h>
#include <LinearFilterBasedSMCAlgorithm.h>
#include <ViterbiAlgorithm.h>
#include <ResamplingCriterion.h>
#include <StdResamplingAlgorithm.h>
#include <StatUtil.h>
#include <Util.h>
#include <lapackpp/gmd.h>
#include <lapackpp/blas1pp.h>
#include <lapackpp/blas2pp.h>
#include <lapackpp/blas3pp.h>
#include <lapackpp/laslv.h>
#include <lapackpp/lavli.h>
#include <mylapack.h>
#include <Particle.h>
#include <ParticleWithChannelEstimation.h>
#include <fstream>
#include <cstdlib>
#include <UnknownChannelOrderAlgorithm.h>

using namespace std;
// using namespace la;

// void Matrix2File(tMatrix A,string name,ofstream f);

tMatrix HsToStackedH(vector<tMatrix> matrices,int m)
{
	
	if((matrices[0].cols() % m)!=0)
		throw RuntimeException("m parameter (memory of the channel) is wrong.");
	int L = matrices[0].rows();
	int N = matrices[0].cols()/m;
	int d = matrices.size()-1;
	
	tMatrix res(matrices[0].rows()*(d+1),N*(m+d));
    res = 0.0;

	for(int i=0;i<=d;i++)
	{
		tRange rowsRange(i*L,(i+1)*L-1);
		tRange colsRange(i*N,i*N+N*m-1);
		res(rowsRange,colsRange).inject(matrices[i]);
	}

	return res;
}

/**
 * 
 * @param argc 
 * @param argv[] 
 * @return 
 */
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

	// --------------------- PAM4 ---------------------------
//     vector<vector<tBit> > secuenciasBits(4,vector<tBit>(2));
//     secuenciasBits[0][0] = 0; secuenciasBits[0][1] = 0;
//     secuenciasBits[1][0] = 0; secuenciasBits[1][1] = 1;
//     secuenciasBits[2][0] = 1; secuenciasBits[2][1] = 0;
//     secuenciasBits[3][0] = 1; secuenciasBits[3][1] = 1;
//     vector<tSymbol> simbolos(4);
//     simbolos[0] = -3; simbolos[1] = -1; simbolos[2] = 1; simbolos[3] = 3;
//     Alphabet pam4(2,4,secuenciasBits,simbolos);
// 
//     vector<tBit> secuenciaAbuscar(2);
//     secuenciaAbuscar[0] = 1; secuenciaAbuscar[1] = 1;
//     tSymbol simboloDevuelto = pam4[secuenciaAbuscar];
	// ----------------------------------------------------------

	// -------------- PAM2 -------------------
//     vector<vector<tBit> > secuenciasBits(2,vector<tBit>(1));
// 	secuenciasBits = *(new vector<vector<tBit> >(2,vector<tBit>(1)));
	vector<vector<tBit> > secuenciasBits(2,vector<tBit>(1));
    secuenciasBits[0][0] = 0; secuenciasBits[1][0] = 1;
//     vector<tSymbol> simbolos(2);
// 	simbolos = *(new vector<tSymbol>(2));
	vector<tSymbol> simbolos(2);
    simbolos[0] = -1; simbolos[1] = 1;
    Alphabet pam2(1,2,secuenciasBits,simbolos);


	vector<tBit> secuenciaDevuelta;
// 	secuenciaDevuelta = pam4[-1];
// 	cout << simboloDevuelto << endl;
// 	cout << "el tamaño de la secuencia es " << secuenciaDevuelta[0];
// 	secuenciaDevuelta.resize(6);
// 	vector<tSymbol> secuenciaSimbolos(8);
// 	pam4.IntToSymbolsArray(13,secuenciaSimbolos);
// 	cout << "Secuencia devuelta" << endl;
// 	for(int i=0;i<secuenciaSimbolos.size();i++)
// 		cout << secuenciaSimbolos[i];
// 	cout << endl;
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
	int L=3,N=2,m=2,K=10;
	int longSecEntr = 5;
	int nParticles = 30;
	int d = m -1;
	double forgettingFactor = 0.9;
	double channelMean=0.0,channelVariance=1.0,ARvariance=0.0001;
	double samplingVariance = 0.0625;
	vector<double> ARcoefficients(1);
	ARcoefficients[0] = 0.99999;
	tRange todasFilasSimbolos(0,N-1);
	double muLMS = 0.05;

	Bits bitsTransmitir(N,K);

	tMatrix simbolosTransmitir = Modulator::Modulate(bitsTransmitir,pam2);

	tMatrix preambulo(N,m-1);
	preambulo = -1.0;
	simbolosTransmitir = Util::Append(preambulo,simbolosTransmitir);

	ARchannel canal(N,L,m,simbolosTransmitir.cols(),channelMean,channelVariance,ARcoefficients,ARvariance);
	cout << "longitud canal " << canal.Length() << " simbolos.co " << simbolosTransmitir.cols() << endl;
	for(int i=canal.Memory()-1;i<canal.Length();i++)
		cout << canal[i] << endl << "****************" << endl;

	// ---------- modelo apilado -----------
	vector<tMatrix> conjuntoMatrices(d+1);
	for(int i=0;i<=d;i++)
	{
		conjuntoMatrices[i] = canal[canal.Memory()-1+i];
		cout << "++++++++++++" << endl << conjuntoMatrices[i] << endl << "++++++++++++++" << endl ;
	}

	tMatrix matrizApilada = HsToStackedH(conjuntoMatrices,m);
	cout << "La matriz apilada es" << endl << matrizApilada << endl;
	// -------------------------------------

	ChannelDependentNoise ruido(&canal);
	ruido.SetSNR(12,1);

	ofstream fid("sal", ofstream::out);
	ofstream fid2("sal2", ofstream::out);

	Util::MatrixToStream(preambulo,"preambulo",fid);
	Util::MatricesVectorToStream(conjuntoMatrices,"matrices",fid2);

	fid.close();
	fid2.close();

// 	for(int iVarianza=31;iVarianza<300;iVarianza++)
// 	{
// 	cout << "varianza = " << ruido.VarianceAt(iVarianza) << " ruido en " << iVarianza << endl << ruido[iVarianza] << endl;
// 	}
// 	exit(1);

// 	cout << "El ruido" << endl;
// 	ruido.Print();
// 
//  for(int iVarianza=0;iVarianza<ruido.Length();iVarianza++)
//  {
//  cout << ruido.VarianceAt(iVarianza) << " , ";
//  }
// 
//     ruido.SetSNR(3,1);
// 
//     cout << "El ruido" << endl;
//     ruido.Print();
// 
//  for(int iVarianza=0;iVarianza<ruido.Length();iVarianza++)
//  {
//  cout << ruido.VarianceAt(iVarianza) << " , ";
//  }
// 
//     ruido.SetSNR(12,1);
// 
//     cout << "El ruido" << endl;
//     ruido.Print();
// 
//  for(int iVarianza=0;iVarianza<ruido.Length();iVarianza++)
//  {
//  cout << ruido.VarianceAt(iVarianza) << " , ";
//  }

    /*
    ofstream f("salida",ofstream::trunc);
    return 0;*/

	ChannelDependentNoise ruidoCopia = ruido;

	tMatrix observaciones = canal.Transmit(simbolosTransmitir,ruidoCopia);

// 	ARchannel canalCopia = canal;

// 	tMatrix observacionesCopia = canalCopia.Transmit(simbolosTransmitir,ruido);

// 	cout << "Las observaciones" << endl << observaciones << endl;
// 
// 	cout << "Las observaciones copia" << endl << observacionesCopia << endl;

	tMatrix mediaInicial(L,N*m);
	mediaInicial = 0.0;

	tMatrix mediaInicial3(L,N*3);
	mediaInicial3 = 0.0;

	tMatrix mediaInicial5(L,N*5);
	mediaInicial3 = 0.0;

	KalmanEstimator estimador(mediaInicial,ARcoefficients[0],ARvariance);

	vector<ChannelMatrixEstimator *> estimadores;

	estimadores.push_back(new KalmanEstimator(mediaInicial,ARcoefficients[0],ARvariance));
	estimadores.push_back(new KalmanEstimator(mediaInicial3,ARcoefficients[0],ARvariance));
	estimadores.push_back(new KalmanEstimator(mediaInicial5,ARcoefficients[0],ARvariance));

	UnknownChannelOrderAlgorithm pr("Orden de canal desconocido",pam2,L,N,K,estimadores,tMatrix(N,10),m-1);	

    return 0;
	
	RLSEstimator estimadorRLS(mediaInicial,forgettingFactor);

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

		RLSEstimator* copiaEstimadorRLS = new RLSEstimator(estimadorRLS);

    RMMSEDetector rmmseDetector(L*(d+1),N*(m+d),1.0,0.98,N*(d+1));

	for(int i=m-1;i<observaciones.cols()-d;i++)
	{
// 		cout << "OOOOOOOOOOOOO i es " << i << " 000000000000000000" << endl;
		tRange rangoColumnas(i-m+1,i);
		tMatrix subMatrizSimbolos = simbolosTransmitir(todasFilasSimbolos,rangoColumnas);
		tMatrix est = estimadorRLS.NextMatrix(observaciones.col(i),subMatrizSimbolos,ruido.VarianceAt(i));
// 		cout << "Estimacion (varianza es " << ruido.VarianceAt(i) << ")" << endl << est << endl;
// 		cout << "El ruido" << endl << ruido
// 		cout << "Canal de verdad" << endl << canal[i] << endl << "-------------" << endl;

// 		RLSEstimator* ptrEstimadorRLS = new RLSEstimator(estimadorRLS);

// 		ParticleWithChannelEstimation part(0.1,10,2,3,4,ptrEstimadorRLS);
   
		tRange rangoSuavizado(i,i+d);
		tRange filasObservaciones(0,L-1);
		tVector observacionesApiladas = Util::ToVector(observaciones(filasObservaciones,rangoSuavizado),columnwise);
// 		cout << "El vector de observacioes apilado concatenado" << endl << observacionesApiladas << endl << "oooooooooooooo" << endl;
//         cout << "correspondiente a las observaciones" << endl;
        for(int ii=i;ii<=i+d;ii++)
            cout << observaciones.col(ii) << endl;
		tMatrix canalApilado = HsToStackedH(canal.Range(i,i+d),m);
//         cout << "La matriz de canal apilada es" << endl << canalApilado << endl;
//         cout << "correspondiente a las matrices" << endl;
//         for(int ii=i;ii<=i+d;ii++)
//             cout << canal[ii] << endl;        


        tVector estimacionesBlandas = rmmseDetector.Detect(observacionesApiladas,canalApilado);
//         cout << "Las estimacionees blandas" << endl << estimacionesBlandas << endl;

        tMatrix matrizRuidoApilado(L,d+1);
        for(int aux = i;aux<=i+d;aux++)
        {
            tVector auxRuido = ruido[aux];
            for(int aux2=0;aux2<L;aux2++)
                matrizRuidoApilado(aux2,aux-i) = auxRuido(aux2);
        }
        tVector ruidoApilado = Util::ToVector(matrizRuidoApilado,columnwise);
//         cout << "El ruido apilado" << endl << ruidoApilado << endl;
//         cout << "correspondiente a" << endl;
//         for(int aux=i;aux<=i+d;aux++)
//             cout << ruido[aux] << endl;

        tRange rangoSimbolosObservacionApilada(i-m+1,i+d);
        tVector simbolosObservacionApilada = Util::ToVector(simbolosTransmitir(todasFilasSimbolos,rangoSimbolosObservacionApilada),columnwise);
//         cout << "Los simbolos apilados" << endl << simbolosObservacionApilada << endl;
//         cout << "correspondiente a" << endl;
//         for(int aux=i-m+1;aux<=i+d;aux++)
//             cout << simbolosTransmitir.col(aux) << endl;

        tVector observacionApiladaConstruida(L*(d+1));
        Blas_Mat_Vec_Mult(canalApilado,simbolosObservacionApilada,observacionApiladaConstruida);
        Util::Add(observacionApiladaConstruida,ruidoApilado,observacionApiladaConstruida);

//         cout << "El vector de observaciones construido" << endl << observacionApiladaConstruida << endl;

// 		vector<tMatrix> conjMatrices = canal.Range(i,i+d);
// 		for(int auxI=0;auxI<conjMatrices.size();auxI++)
// 			cout << conjMatrices[auxI] << "d es " << d << " auxI es " << auxI << " numero de matrices " << conjMatrices.size() << endl;

/*		double verosimil = estimador.Likelihood(observaciones.col(i),subMatrizSimbolos,ruido.VarianceAt(i));
		cout << "La verosimilitud=" << verosimil << endl;

		tMatrix simbolosChungos(N,m);
		simbolosChungos = -11;
		double verosimilChunga = estimador.Likelihood(observaciones.col(i),simbolosChungos,ruido.VarianceAt(i));
		cout << "La verosimilitud chunga=" << verosimilChunga << endl;	*/	
	}

// 	tVector hola1(2);
// 	hola1(0) = 1.2; hola1(1) = 0.3;
// 
// 	tVector hola2(4);
// 	hola2(0) = 11.2; hola2(1) = 3.1; hola2(2) = 0.1; hola2(3) = 0.1001;
// 
// 	tMatrix resHola(2,4);
// 	
// 	Util::Mult(hola1,hola2,resHola);
// 
// 	cout << "hola1" << endl << hola1 << endl << "hola2" << endl << hola2 << endl;
// 	cout << "res" << endl << resHola << endl;
	

// ML_SMCAlgorithm::ML_SMCAlgorithm(string name, Alphabet alphabet,int L,int N, ChannelMatrixEstimator& channelEstimator, tMatrix preamble, int smoothingLag, int nParticles, ResamplingCriterion resamplingCriterion): SMCAlgorithm(name, alphabet, channelEstimator, preamble, smoothingLag, nParticles, resamplingCriterion)

// 	tMatrix preambulo(N,m-1);
// 	preambulo = -1.0;

	ResamplingCriterion criterioRemuestreo(0.9);
	StdResamplingAlgorithm algoritmoRemuestreo;

	ML_SMCAlgorithm algoritmo("Detector suavizado optimo",pam2,L,N,K-d,&estimador,preambulo,m-1,nParticles,criterioRemuestreo,algoritmoRemuestreo);

// 	cout << "El canal en pruebas" << endl << canal[55] << endl;
// 
// cout << "El canal en pruebas" << endl << canal[55] << endl;

	RMMSEDetector detectorMMSE(L*(d+1),N*(m+d),1.0,0.98,N*(d+1));

	RLSEstimator estimadorRLSfiltroLineal(mediaInicial,forgettingFactor);
	LMSEstimator estimadorLMSfiltroLineal(mediaInicial,muLMS);

	LinearFilterBasedSMCAlgorithm algoritmoFiltroLineal("Filtro lineal",pam2,L,N,K-d,&estimadorRLSfiltroLineal,&detectorMMSE,preambulo,m-1,nParticles,criterioRemuestreo,algoritmoRemuestreo,ARcoefficients[0],samplingVariance,ARvariance);

// 	cout << "El canal en pruebas" << endl << canal[55] << endl;

// 	char c; cin >> c;


	tMatrix secEntrenamiento = simbolosTransmitir(todasFilasSimbolos,*(new tRange(m-1,m+longSecEntr-2)));
// 	algoritmo.Run(observaciones,ruido.Variances());

// 	algoritmo.Run(observaciones,ruido.Variances(),secEntrenamiento);
// 	double pe = algoritmo.SER(simbolosTransmitir(todasFilasSimbolos,*(new tRange(m-1+longSecEntr,simbolosTransmitir.cols()-d-1))));
// 	cout << "La probabilidad de error es " << pe << endl;
// 	ojo: los ultimos simbolos no se detectan

	algoritmoFiltroLineal.Run(observaciones,ruido.Variances(),secEntrenamiento);
	double pe = algoritmoFiltroLineal.SER(simbolosTransmitir(todasFilasSimbolos,*(new tRange(m-1+longSecEntr,simbolosTransmitir.cols()-d-1))));
	cout << "La probabilidad de error es " << pe << " y el MSE " << algoritmoFiltroLineal.MSE(canal.Range(K*9/10,simbolosTransmitir.cols()-d-1)) << endl;

    cout << "La probabilidad de error es " << pe << " y el MSE " << algoritmoFiltroLineal.MSE(canal.Range(m-1+longSecEntr,simbolosTransmitir.cols()-d-1)) << endl;

//     vector<tMatrix> matricesDetectadas = algoritmoFiltroLineal.GetEstimatedChannelMatrices();
//     for(int i2=0;i2<matricesDetectadas.size();i2++)
//         cout << "Instante " << i2 << " matriz estimada: " << endl << matricesDetectadas[i2] << endl;
// 
//     cout << "Canal de verdad" << endl << canal[observaciones.cols()-d] << endl;
//     cout << "La probabilidad de error es (con SER2) " << algoritmoFiltroLineal.SER2(simbolosTransmitir(todasFilasSimbolos,*(new tRange(m-1+longSecEntr,simbolosTransmitir.cols()-d-1)))) << endl;

//     cout << "Los simbolos detectados son" << endl << algoritmoFiltroLineal.GetDetectedSymbolVectors() << endl;
//     cout << "Los voy a comparar con" << endl << simbolosTransmitir(todasFilasSimbolos,*(new tRange(m-1+longSecEntr,simbolosTransmitir.cols()-d-1))) << endl;
	// ojo: los ultimos simbolos no se detectan

// 	cout << "ahi va" << algoritmo._estimatedChannelMatrices[0][0] << endl;


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

//     tVector v1(4),v2(4);
//     tMatrix C(4,4);
//     v1 = 1.2;
//     v2 = 2.3;
//     Util::Mult(v1,v2,C);
// 
//     cout << "Res es" << endl << C << endl;

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

	// -------------- ordenacion ------------------------
// 	vector<int> indexes(4);
// 	indexes[0] = 4;indexes[1] = 2;indexes[2] = 9;indexes[3] = -1;
// 	for(int holita=0;holita<indexes.size();holita++)
// 		cout << indexes[holita] << "," ;
// 	cout << endl;
// 	vector<int>::iterator startIterator = indexes.begin();
// 	vector<int>::iterator endIterator = indexes.end();
// 	sort(startIterator,endIterator);
// 	for(int holita=0;holita<indexes.size();holita++)
// 		cout << indexes[holita] << "," ;
// 	cout << endl;
	// ----------------------------------------------------

	// --------------------- Remuestreo -------------------------------
// 	tMatrix **estimatedChannelMatrices;
// 	tMatrix **detectedSymbols;
// 	ChannelMatrixEstimator **particlesChannelMatrixEstimators;
// 	
// 	estimatedChannelMatrices = new tMatrix*[3];
// 	detectedSymbols = new tMatrix*[3];
// 	particlesChannelMatrixEstimators = new ChannelMatrixEstimator*[3];
// 	for(int i=0;i<3;i++)
// 	{
// 		estimatedChannelMatrices[i] = new tMatrix[5];
// 		particlesChannelMatrixEstimators[i] = new KalmanEstimator(ARcoefficients[0],ARvariance,*(new tMatrix(i+1,4)));
// 		for(int j=0;j<5;j++)
// 		{
// 			estimatedChannelMatrices[i][j] = *(new tMatrix(2,2));
// 		}
// 		detectedSymbols[i] = new tMatrix(2,4);
// 		(*detectedSymbols[i])(0,0) = i;
// 		estimatedChannelMatrices[i][0](0,0) = i;
// 	}
// 
// 	estimatedChannelMatrices[0][0](0,0) = 1.1;
// 
// 	//impresion
// 	for(int i=0;i<3;i++)
// 	{
// 		cout << "Matrices de canal estimadas" << endl;
// 		for(int j=0;j<5;j++)
// 		{
// 			cout << estimatedChannelMatrices[i][j];
// 		}
// 		cout << "Simbolos" << endl << *detectedSymbols[i] << endl;
// 		cout << "Filas del estimador=" << (particlesChannelMatrixEstimators[i])->Rows() << endl;
// 		cout << "---------" << endl;
// 	}
// 
// 	vector<int> indices(3);
// 	indices[0] = 1;indices[1] = 0;indices[2] = 2;
// 
// 	StdResamplingAlgorithm::Resampling(&estimatedChannelMatrices,&detectedSymbols,&particlesChannelMatrixEstimators,indices,3,0,4,5);
// 
// 	cout << "Despues del remuestreo..." << endl << endl;
// 
// 	//impresion
// 	for(int i=0;i<3;i++)
// 	{
// 		cout << "Matrices de canal estimadas" << endl;
// 		for(int j=0;j<5;j++)
// 		{
// 			cout << estimatedChannelMatrices[i][j];
// 		}
// 		cout << "Simbolos" << endl << *detectedSymbols[i] << endl;
// 		cout << "Filas del estimador=" << (particlesChannelMatrixEstimators[i])->Rows() << endl;
// 		cout << "---------" << endl;
// 	}
	//----------------------------------------------------------------------


	// ---------------------- Maximo --------------------------------
// 	tVector pruebaMaximo(4);
// 	pruebaMaximo(0) = -10.2;pruebaMaximo(1)=2.132;pruebaMaximo(2)=234.1;pruebaMaximo(3) = -12.1;
// 	int imax;
// 	Util::Max(pruebaMaximo,imax);
// 	cout << "El vector" << endl << pruebaMaximo << endl;
// 	cout << "max=" << pruebaMaximo(imax) << " indice=" << imax << endl;
	// --------------------------------------------------------------

	// --------------- RandnArray ----------------------------
// 	Random r;
// 	cout << "La matriz aleatoria" << endl << StatUtil::RandnMatrix(2,3,0.0,1.0,r) << endl;
	// --------------------------------------------------------

// 	cout << "1º simbolo: " << pam2[0] << " 2º: " << pam2[1] << endl;
// 	cout << "Normal= " << StatUtil::NormalPdf(3.1,2.0,1.0) << endl;
// 	int tam = 8;
// 	tVector x(tam);
// 	x = 2.0;
// 	tVector mean(tam);
// 	mean = 1.0;
// 	tMatrix covarianza = StatUtil::RandnMatrix(tam,tam,0.0,1.0);
// 	cout << "la covarianza es" << endl << covarianza << endl;
// 	cout << "El resultado=" << StatUtil::NormalPdf(x,mean,covarianza);

// 	tMatrix holita = StatUtil::RandnMatrix(2,4,0.0,1.0);
// 	cout << "holita es" << endl << holita << endl;
// 	tVector holitaVector = Util::ToVector(holita,columnwise);
// 	cout << "como vector" << endl << holitaVector << endl;
// 	cout << "vuelve a matriz" << endl << Util::ToMatrix(holitaVector,columnwise,2,8) << endl;

	// ********************* Pruebas con particulas ****************************
// 	Particle particula(0.3,3,5);
// 	Particle particula2(0.3,3,5);
// 
// 	particula2 = particula;
// 
// 	tVector v(3);
// 	v(0) = 1.1; v(1) = 2.123413; v(2) = 3.00001;
// 
// 	particula.Print();
// 
// 	particula.SetWeight(0.314);
// 
// 	particula.SetSymbolVector(2,v);
// 
// 	tMatrix pruebita = particula.GetSymbolVectors(1,3);
// 
// 	particula.Print();
// 
// 	cout << "un rango" << endl << pruebita << endl;
// 
// 	ParticleWithChannelEstimation part(0.1,2,10,3,4,copiaEstimadorRLS);
// 
// 	tMatrix matrizInsertar = StatUtil::RandnMatrix(3,4,0.0,1.0);
// 
// 	cout << "Matriz que se va a insertar" << endl << matrizInsertar << endl;
// 
// 	part.SetChannelMatrix(6,matrizInsertar);
// 
// 	cout << part.GetChannelMatrix(6) << endl;
// 
// 	ParticleWithChannelEstimation part2 = part;
// 	cout << part2.GetChannelMatrix(6) << endl;
// 
// 	part2.SetChannelMatrix(3,matrizInsertar);
// 
// 	part = part2;
// 
// 	cout << "despues del igual" << endl << part.GetChannelMatrix(3) << endl;
	// ***********************************************************************

// 	int a** = new int[][4];
// 	int a = 4;
// 	vector< vector<int> > hola(a);

	cout << "Al final del programa" << endl << endl;

	ViterbiAlgorithm algoritmoViterbi("Viterbi",pam2,L,N,K-d,canal,preambulo,d);
	algoritmoViterbi.Run(observaciones,ruido.Variances());
// 	algoritmoViterbi.PrintStage(exitStage);
// 	cout << "Prob error es " << algoritmoViterbi.SER(simbolosTransmitir(todasFilasSimbolos,*(new tRange(m-1+longSecEntr,simbolosTransmitir.cols()-1)))) << endl;

    cout << "Prob error es " << algoritmoViterbi.SER(simbolosTransmitir(todasFilasSimbolos,*(new tRange(m-1+longSecEntr,simbolosTransmitir.cols()-d-1)))) << endl;

//     cout << "Prob error es (con SER2) " << algoritmoViterbi.SER2(simbolosTransmitir(todasFilasSimbolos,*(new tRange(m-1+longSecEntr,simbolosTransmitir.cols()-d-1)))) << endl;

//     cout << "Los simbolos detectados son" << endl << algoritmoViterbi.GetDetectedSymbolVectors() << endl;
//     cout << "Los voy a comparar con" << endl << simbolosTransmitir(todasFilasSimbolos,*(new tRange(m-1+longSecEntr,simbolosTransmitir.cols()-d-1))) << endl;

// m-1+longSecEntr,simbolosTransmitir.cols()-d-1

	// ********************** SymbolsArrayToInt **********************
// 	vector<tSymbol> vectorSimb(4,1);
// 	pam2.IntToSymbolsArray(7,vectorSimb);
// 	for(int holita=0;holita<4;holita++)
// 		cout << vectorSimb[holita] << " ";
// 	cout << endl;
// 	cout << "Covertido a numero " << pam2.SymbolsArrayToInt(vectorSimb) << endl;
	// **********************************************************

// 	cout << simbolosTransmitir << endl;

//     vector<tMatrix> holita;
//     cout << "El tamaño de holita es " << holita.size() << endl;
// 
//     tMatrix prueba1(2,2);
//     prueba1(0,0) = 1.0;prueba1(0,1) = 2.0;prueba1(1,0) = 3.0;prueba1(1,1) = 4.0;
// 
//     tMatrix prueba2 = prueba1;
//     prueba2(0,0) = 3.0;
// 
//     cout << Util::SquareError(prueba1,prueba2) << endl;

    return 0;
}
