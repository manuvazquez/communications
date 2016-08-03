#ifndef KNOWNSYMBOLSNOBSERVATIONSKALMANESTIMATOR_H
#define KNOWNSYMBOLSNOBSERVATIONSKALMANESTIMATOR_H

#include <KalmanEstimator.h>

/**
    @author Manu <manu@rustneversleeps>
*/
class KnownSymbolsNObservationsKalmanEstimator : public KalmanEstimator
{
protected:
    uint _presentTime;
    const MatrixXd &_symbols;
	const MatrixXd &_observations;
public:
    KnownSymbolsNObservationsKalmanEstimator(const MatrixXd& initialEstimation, const MatrixXd& variances, uint N, vector< double > ARcoefficients, double ARvariance, const MatrixXd& symbols, const MatrixXd& observations, uint startDetectionTime);

    KnownSymbolsNObservationsKalmanEstimator* clone() const;

    virtual MatrixXd nextMatrix(const VectorXd &observations, const MatrixXd &symbolsMatrix, double noiseVariance);


};

#endif
