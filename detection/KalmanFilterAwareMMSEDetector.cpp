/*
    Copyright 2012 Manu <manuavazquez@gmail.com>

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
*/


#include "KalmanFilterAwareMMSEDetector.h"

#include <math.h>
#include <TransmissionUtil.h>

// #define DEBUG
// #define DEBUG2
// #define DEBUG3

KalmanFilterAwareMMSEDetector::KalmanFilterAwareMMSEDetector(uint rows, uint cols, double alphabetVariance,uint nSymbolsToBeDetected,KalmanEstimator *kalmanEstimator,std::vector<double> ARcoefficients)
:MMSEDetector(rows,cols,alphabetVariance,nSymbolsToBeDetected),_kalmanEstimator(kalmanEstimator),_ARcoefficients(ARcoefficients)
{
	// this implementation is only meaningful for channels that evolve according with a first order AR process
	assert(ARcoefficients.size()==1);
}

VectorXd KalmanFilterAwareMMSEDetector::detect(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance)
{
	uint nRows = channelMatrix.rows();
	
#ifdef DEBUG
	cout << "channelMatrix.cols() = " << channelMatrix.cols() << " nRows =  " << nRows << endl;
	cout << "channelMatrix:" << endl << channelMatrix << endl;
	cout << "_kalmanEstimator->predictedMatrix(): " << endl << _kalmanEstimator->predictedMatrix() << endl;
#endif
	
	uint N = _kalmanEstimator->nInputs();
	int Nm = _kalmanEstimator->cols();
	uint m = uint(Nm)/N;
	uint L = _kalmanEstimator->rows();
	
	// smoothing lag "d" is inferred from the received channel matrix and the internal Kalman filter
	uint d = nRows/_kalmanEstimator->rows() -1;
	
#ifdef DEBUG
	cout << "m = " << _kalmanEstimator->memory() << " d = " << d << " Nm = " << Nm << endl;
	cout << "_kalmanEstimator->getPredictiveCovariance()" << endl << _kalmanEstimator->getPredictiveCovariance() << endl;
#endif
	
	std::vector<MatrixXd> predictedMatrices(d+1);
	std::vector<MatrixXd> predictedCovariances(d+1);
	
	predictedMatrices[0] = _kalmanEstimator->predictedMatrix();
	predictedCovariances[0] = _kalmanEstimator->getPredictiveCovariance();
	
	KalmanEstimator *kalmanEstimatorClone = _kalmanEstimator->clone();
	
	for(uint i=1;i<=d;i++)
	{
		// one PREDICTIVE step is advanced on the cloned Kalman filter: observations vector and symbols matrix are zero so that the filtered mean and covariance be equal to the predictive ones.
		// Regarding the noise variance, a value different from 0 (any one should do) must be passed or otherwise, the matrix inversion within the KF will give rise to NaN's
		kalmanEstimatorClone->nextMatrix(VectorXd::Zero(L),MatrixXd::Zero(N,m),1.0);
		
		predictedMatrices[i] = kalmanEstimatorClone->predictedMatrix();
		predictedCovariances[i] = kalmanEstimatorClone->getPredictiveCovariance();

#ifdef DEBUG
		cout << "predicted for i = " << i << endl << predictedMatrices[i] << endl;
		cout << "covariance: " << endl << predictedCovariances[i] << endl;
#endif
	}
	
	delete kalmanEstimatorClone;
	
	MatrixXd predictedStackedChannelMatrix = TransmissionUtil::channelMatrices2stackedChannelMatrix(predictedMatrices,m);
	
#ifdef DEBUG
	cout << "predictedStackedChannelMatrix:" << endl << predictedStackedChannelMatrix << endl;
#endif
	
	MatrixXd columnsAutoCorrelationSum = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
	
#ifdef DEBUG2
	MatrixXd onlyCovariances = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
#endif
	
	
	// for every possible column index, the indexes within the KF state vector that give that column are obtained
	std::vector<std::vector<uint> > iCol2indexesWithinKFstateVector(Nm);
	
	for(int iCol=0;iCol<Nm;iCol++)
	{
		iCol2indexesWithinKFstateVector[iCol] = _kalmanEstimator->colIndexToIndexesWithinKFstateVector(iCol);
#ifdef DEBUG
		cout << "iCol2indexesWithinKFstateVector for iCol =  " << iCol << endl << iCol2indexesWithinKFstateVector[iCol] << endl;
#endif
	}
	
	for(uint iCol=0;iCol<predictedStackedChannelMatrix.cols();iCol++)
	{
#ifdef DEBUG
		cout << "============= iCol = " << iCol << " ==============" << endl;
#endif
		
		MatrixXd thisColumnCovariance = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
		
		for(uint i=0;i<(d+1);i++)
		{
			int upperSubcolumnOriginalColumn = iCol - (i*N);
			for(uint j=i;j<(d+1);j++)
			{
				int lowerSubcolumnOriginalColumn = iCol - (j*N);
				
				if(upperSubcolumnOriginalColumn>=Nm || upperSubcolumnOriginalColumn<0)
					continue;
				
				if(lowerSubcolumnOriginalColumn>=Nm || lowerSubcolumnOriginalColumn<0)
					continue;
				
#ifdef DEBUG
				cout << "upperSubcolumnOriginalColumn = " << upperSubcolumnOriginalColumn << " lowerSubcolumnOriginalColumn = " << lowerSubcolumnOriginalColumn << endl;
#endif
				
				if(i == j) // => upperSubcolumnOriginalColumn == lowerSubcolumnOriginalColumn
					thisColumnCovariance.block(i*L,j*L,L,L) = Util::subMatrixFromVectorIndexes(predictedCovariances[j],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn]);
				else
				{
					MatrixXd subCovariance = Util::subMatrixFromVectorIndexes(predictedCovariances[0],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn])
												*pow(_ARcoefficients[0],double(j-i));
					thisColumnCovariance.block(i*L,j*L,L,L) = subCovariance;
					thisColumnCovariance.block(j*L,i*L,L,L) = subCovariance.transpose(); // ...due to the symmetry of the covariance matrix
				}
			} // for(uint j=i;j<(d+1);j++)
		}
		
		columnsAutoCorrelationSum += thisColumnCovariance + predictedStackedChannelMatrix.col(iCol)*predictedStackedChannelMatrix.col(iCol).transpose();
// 		columnsAutoCorrelationSum += predictedStackedChannelMatrix.col(iCol)*predictedStackedChannelMatrix.col(iCol).transpose();
		
#ifdef DEBUG2
		onlyCovariances += thisColumnCovariance;
#endif
		
#ifdef DEBUG3
		cout << "thisColumnCovariance" << endl << thisColumnCovariance << endl;
		cout << "predictedCovariances[0]" << endl << predictedCovariances[0] << endl;
		getchar();
#endif
		
#ifdef DEBUG
		cout << "thisColumnCovariance" << endl << thisColumnCovariance << endl;
		getchar();
#endif
	} // for(uint iCol=0;iCol<predictedStackedChannelMatrix.cols();iCol++)
	
#ifdef DEBUG2
	cout << "onlyCovariances" << endl << onlyCovariances << endl;
	getchar();
#endif

	MatrixXd _Rx = noiseCovariance + _alphabetVariance*columnsAutoCorrelationSum;

    _filter = _Rx.inverse()*channelMatrix*_alphabetVariance;

    VectorXd softEstimations = _filter.transpose()*observations;

    // required for nthSymbolVariance computing
    _channelMatrix = channelMatrix;

    return softEstimations.segment(_detectionStart,_nSymbolsToBeDetected);
}

KalmanFilterAwareMMSEDetector* KalmanFilterAwareMMSEDetector::clone()
{
	return new KalmanFilterAwareMMSEDetector(*this);
}
