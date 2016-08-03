#include "KnownSymbolsNObservationsKalmanEstimator.h"

// #define DEBUG

KnownSymbolsNObservationsKalmanEstimator::KnownSymbolsNObservationsKalmanEstimator(const MatrixXd& initialEstimation, const MatrixXd& variances, uint N, vector<double> ARcoefficients, double ARvariance,const MatrixXd &symbols,const MatrixXd &observations,uint startDetectionTime): KalmanEstimator(initialEstimation, variances, N, ARcoefficients, ARvariance),_presentTime(startDetectionTime),_symbols(symbols),_observations(observations)
{
}

MatrixXd KnownSymbolsNObservationsKalmanEstimator::nextMatrix(const VectorXd &observations, const MatrixXd &symbolsMatrix, double noiseVariance)
{
    _presentTime++;
	
#ifdef DEBUG
	std::cout << "_observations.col(_presentTime-1):" << std::endl << _observations.col(_presentTime-1) << std::endl;
#endif
	
	return KalmanEstimator::nextMatrix(_observations.col(_presentTime-1), _symbols.block(0,_presentTime-_channelOrder,_nInputs,_channelOrder), noiseVariance);
}

KnownSymbolsNObservationsKalmanEstimator* KnownSymbolsNObservationsKalmanEstimator::clone() const
{
	return new KnownSymbolsNObservationsKalmanEstimator(*this);
}
