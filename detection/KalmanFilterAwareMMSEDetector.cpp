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
#include <map>
#include <TransmissionUtil.h>

// #define DEBUG

KalmanFilterAwareMMSEDetector::KalmanFilterAwareMMSEDetector(uint rows, uint cols, double alphabetVariance,uint nSymbolsToBeDetected,KalmanEstimator *kalmanEstimator,std::vector<double> ARcoefficients,bool interferenceCancellation)
:MMSEDetector(rows,cols,alphabetVariance,nSymbolsToBeDetected),_kalmanEstimator(kalmanEstimator),_ARcoefficients(ARcoefficients),_interferenceCancellation(interferenceCancellation)
{
}

VectorXd KalmanFilterAwareMMSEDetector::detect(const VectorXd &observations, const MatrixXd &channelMatrix, const MatrixXd& noiseCovariance)
{	
	if(_ARcoefficients.size()>1)
		return detect2orderAndAboveARprocess(observations,channelMatrix,noiseCovariance);
	
	uint nRows = channelMatrix.rows();
	
	uint N = _kalmanEstimator->nInputs();
	int Nm = _kalmanEstimator->cols();
	uint m = uint(Nm)/N;
	uint L = _kalmanEstimator->rows();
	
	// index of the first stacked channel matrix column to be accounted for when computing the channel matrix covariance
	uint iFirstColumn;
	
	if(_interferenceCancellation)
		iFirstColumn = N*(m-1);
	else
		iFirstColumn = 0;
	
	// smoothing lag "d" is inferred from the received channel matrix and the internal Kalman filter
	uint d = nRows/_kalmanEstimator->rows() -1;
	
	std::vector<MatrixXd> predictedMatrices(d+1);
	std::vector<MatrixXd> predictedCovariances(d+1);
	
	predictedMatrices[0] = _kalmanEstimator->predictedMatrix();
	predictedCovariances[0] = _kalmanEstimator->getInternalPredictiveCovariance();

	KalmanEstimator *kalmanEstimatorClone = _kalmanEstimator->clone();
	
	for(uint i=1;i<=d;i++)
	{
		// one PREDICTIVE step is advanced on the cloned Kalman filter: observations vector and symbols matrix are zero so that the filtered mean and covariance be equal to the predictive ones.
		// As for the noise variance, a value different from 0 (any should do) must be passed or otherwise, the matrix inversion within the KF will give rise to NaN's
		kalmanEstimatorClone->nextMatrix(VectorXd::Zero(L),MatrixXd::Zero(N,m),1.0);
		
		predictedMatrices[i] = kalmanEstimatorClone->predictedMatrix();
		predictedCovariances[i] = kalmanEstimatorClone->getInternalPredictiveCovariance();
	}
	
	delete kalmanEstimatorClone;
	
	MatrixXd predictedStackedChannelMatrix = TransmissionUtil::channelMatrices2stackedChannelMatrix(predictedMatrices,m);

	MatrixXd columnsAutoCorrelationSum = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());

	// for every possible column index, the indexes within the KF state vector that give that column are obtained
	std::vector<std::vector<uint> > iCol2indexesWithinKFstateVector(Nm);
	
	for(int iCol=0;iCol<Nm;iCol++)
		iCol2indexesWithinKFstateVector[iCol] = _kalmanEstimator->colIndexToIndexesWithinKFstateVector(iCol);
	
	for(uint iCol=iFirstColumn;iCol<predictedStackedChannelMatrix.cols();iCol++)
	{
		MatrixXd thisColumnCovariance = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
		
		for(uint iUpperSubcolumn=0;iUpperSubcolumn<(d+1);iUpperSubcolumn++)
		{
			int upperSubcolumnIndexWithinItsMatrix = iCol - (iUpperSubcolumn*N);
			for(uint iLowerSubcolumn=iUpperSubcolumn;iLowerSubcolumn<(d+1);iLowerSubcolumn++)
			{
				int lowerSubcolumnIndexWithinItsMatrix = iCol - (iLowerSubcolumn*N);
				
				// correlation with a vector of zeros is zero
				if(upperSubcolumnIndexWithinItsMatrix>=Nm || upperSubcolumnIndexWithinItsMatrix<0)
					continue;
				
				if(lowerSubcolumnIndexWithinItsMatrix>=Nm || lowerSubcolumnIndexWithinItsMatrix<0)
					continue;
				
				if(iUpperSubcolumn == iLowerSubcolumn) // => upperSubcolumnIndexWithinItsMatrix == lowerSubcolumnIndexWithinItsMatrix
					thisColumnCovariance.block(iUpperSubcolumn*L,iUpperSubcolumn*L,L,L) = Util::subMatrixFromVectorIndexes(predictedCovariances[iUpperSubcolumn],iCol2indexesWithinKFstateVector[upperSubcolumnIndexWithinItsMatrix],iCol2indexesWithinKFstateVector[upperSubcolumnIndexWithinItsMatrix]);
				else
				{
					// only "predictedCovariances[0]" is needed because the difference between a column in a given time instant and the same column a number of time steps later is given by a factor (which is the power below)
					MatrixXd subCovariance = Util::subMatrixFromVectorIndexes(predictedCovariances[0],iCol2indexesWithinKFstateVector[upperSubcolumnIndexWithinItsMatrix],iCol2indexesWithinKFstateVector[lowerSubcolumnIndexWithinItsMatrix])
												*pow(_ARcoefficients[0],double(iUpperSubcolumn+iLowerSubcolumn));
					thisColumnCovariance.block(iUpperSubcolumn*L,iLowerSubcolumn*L,L,L) = subCovariance;
					thisColumnCovariance.block(iLowerSubcolumn*L,iUpperSubcolumn*L,L,L) = subCovariance.transpose(); // ...due to the symmetry of the covariance matrix
				}
			} // for(uint j=i;j<(d+1);j++)
		}
		
		columnsAutoCorrelationSum += thisColumnCovariance + predictedStackedChannelMatrix.col(iCol)*predictedStackedChannelMatrix.col(iCol).transpose();
// 		columnsAutoCorrelationSum += predictedStackedChannelMatrix.col(iCol)*predictedStackedChannelMatrix.col(iCol).transpose();

	} // for(uint iCol=0;iCol<predictedStackedChannelMatrix.cols();iCol++)

	MatrixXd _Rx = noiseCovariance + _alphabetVariance*columnsAutoCorrelationSum;

    _filter = _Rx.inverse()*channelMatrix*_alphabetVariance;

    VectorXd softEstimations = _filter.transpose()*observations;

    // required for nthSymbolVariance computing
    _channelMatrix = channelMatrix;

    return softEstimations.segment(_detectionStart,_nSymbolsToBeDetected);
}

VectorXd KalmanFilterAwareMMSEDetector::detect2orderAndAboveARprocess(const VectorXd &observations, const MatrixXd &channelMatrix, const MatrixXd& noiseCovariance)
{
	uint nRows = channelMatrix.rows();
	
	uint N = _kalmanEstimator->nInputs();
	uint Nm = _kalmanEstimator->cols();
	uint m = uint(Nm)/N;
	uint L = _kalmanEstimator->rows();
	
	// index of the first stacked channel matrix column to be accounted for when computing the channel matrix covariance
	uint iFirstColumn;
	
	if(_interferenceCancellation)
		iFirstColumn = N*(m-1);
	else
		iFirstColumn = 0;
	
	// smoothing lag "d" is inferred from the received channel matrix and the internal Kalman filter
	uint d = nRows/_kalmanEstimator->rows() -1;
	
	// TODO: if care is taken, there is no need to compute "predictedStackedChannelMatrix" because "channelMatrix" should contain all the relevant information
	
	std::vector<MatrixXd> predictedMatrices(d+1);
	std::vector<MatrixXd> predictedCovariances(d+1);
	
	predictedMatrices[0] = _kalmanEstimator->predictedMatrix();
	predictedCovariances[0] = _kalmanEstimator->getInternalPredictiveCovariance();

	KalmanEstimator *kalmanEstimatorClone = _kalmanEstimator->clone();
	
	for(uint i=1;i<=d;i++)
	{
		// one PREDICTIVE step is advanced on the cloned Kalman filter: observations vector and symbols matrix are zero so that the filtered mean and covariance be equal to the predictive ones.
		// As for the noise variance, a value different from 0 (any one should do) must be passed or otherwise, the matrix inversion within the KF will give rise to NaN's
		kalmanEstimatorClone->nextMatrix(VectorXd::Zero(L),MatrixXd::Zero(N,m),1.0);
		
		predictedMatrices[i] = kalmanEstimatorClone->predictedMatrix();
		predictedCovariances[i] = kalmanEstimatorClone->getInternalPredictiveCovariance();
	}
	
	delete kalmanEstimatorClone;
	
	MatrixXd predictedStackedChannelMatrix = TransmissionUtil::channelMatrices2stackedChannelMatrix(predictedMatrices,m);

	MatrixXd columnsAutoCorrelationSum = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());

	// for every possible column index, the indexes within the KF state vector that give that column (for the current time instant) are obtained
	std::vector<std::vector<uint> > iCol2indexesWithinKFstateVector(Nm);
	for(uint iCol=0;iCol<Nm;iCol++)
		iCol2indexesWithinKFstateVector[iCol] = _kalmanEstimator->colIndexToIndexesWithinKFstateVector(iCol,_ARcoefficients.size()-1);
	
	// a map where all the "subcovariances" that will be needed later are stored
	std::map<CovarianceId,MatrixXd> covariancesMap;

	// all the "subcovariances" that are contained in the covariance computed by the KF for the present time, t, are obtained ("t1" and "t2" give the time shifts with respect to t, which is represented by time shift of zero)
	for(int t1=0;t1>-int(_ARcoefficients.size());t1--)
		for(int t2=0;t2>-int(_ARcoefficients.size());t2--)
			for(uint i=0;i<Nm;i++)
				for(uint j=0;j<i;j++)
					covariancesMap[CovarianceId(t1,t2,i,j)]= Util::subMatrixFromVectorIndexes(predictedCovariances[0],
													_kalmanEstimator->colIndexToIndexesWithinKFstateVector(i,_ARcoefficients.size()-1+t1),
													_kalmanEstimator->colIndexToIndexesWithinKFstateVector(j,_ARcoefficients.size()-1+t2));
	
					
	// all the "subcovariances" for every combination of "t1" (value of the time instant corresponding to the upper "subcolumn)  belonging to [-R+1,0] and "t2" (value of the time instant for the lower "subcolumn") belonging to [0,d]
	for(int t1 = -(_ARcoefficients.size()-1);t1<=0;t1++)
		for(int t2=1;t2<=int(d);t2++)
			for(uint i=0;i<Nm;i++)
				for(uint j=0;j<i;j++)
				{
					MatrixXd covariance = MatrixXd::Zero(L,L);
					for(uint r=1;r<=_ARcoefficients.size();r++)
					{
						assert(covariancesMap.find(CovarianceId(t1,t2-r,i,j))!=covariancesMap.end());
						covariance += _ARcoefficients[r-1]*covariancesMap[CovarianceId(t1,t2-r,i,j)];
					}
					covariancesMap[CovarianceId(t1,t2,i,j)] = covariance;
				}

	// all the "subcovariances" for every combination of "t1" belonging to [1,d] and "t2" belonging to [t1+1,d]
	for(int t1=1;t1<=int(d);t1++)
		for(int t2=t1+1;t2<=int(d);t2++)
			for(uint i=0;i<Nm;i++)
				for(uint j=0;j<i;j++)
				{
					MatrixXd covariance = MatrixXd::Zero(L,L);
					for(uint r=1;r<=_ARcoefficients.size();r++)
					{
						assert(covariancesMap.find(CovarianceId(t1-r,t2,i,j))!=covariancesMap.end());
						covariance += _ARcoefficients[r-1]*covariancesMap[CovarianceId(t1-r,t2,i,j)];
					}
					covariancesMap[CovarianceId(t1,t2,i,j)] = covariance;
				}
	
	// -------------
	
	for(uint iCol=iFirstColumn;iCol<predictedStackedChannelMatrix.cols();iCol++)
	{
		MatrixXd thisColumnCovariance = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
		
		for(uint iUpperSubcolumn=0;iUpperSubcolumn<(d+1);iUpperSubcolumn++)
		{
			int upperSubcolumnIndexWithinItsMatrix = iCol - (iUpperSubcolumn*N);
			for(uint iLowerSubcolumn=iUpperSubcolumn;iLowerSubcolumn<(d+1);iLowerSubcolumn++)
			{
				int lowerSubcolumnIndexWithinItsMatrix = iCol - (iLowerSubcolumn*N);
				
				// correlation with a vector of zeros is zero
				if(upperSubcolumnIndexWithinItsMatrix>=int(Nm) || upperSubcolumnIndexWithinItsMatrix<0)
					continue;
				
				if(lowerSubcolumnIndexWithinItsMatrix>=int(Nm) || lowerSubcolumnIndexWithinItsMatrix<0)
					continue;
				
				if(iUpperSubcolumn == iLowerSubcolumn) // => upperSubcolumnIndexWithinItsMatrix == lowerSubcolumnIndexWithinItsMatrix
					thisColumnCovariance.block(iUpperSubcolumn*L,iUpperSubcolumn*L,L,L) = Util::subMatrixFromVectorIndexes(predictedCovariances[iUpperSubcolumn],iCol2indexesWithinKFstateVector[upperSubcolumnIndexWithinItsMatrix],iCol2indexesWithinKFstateVector[upperSubcolumnIndexWithinItsMatrix]);
				else
				{
					MatrixXd subCovariance = covariancesMap[CovarianceId(iUpperSubcolumn,iLowerSubcolumn,upperSubcolumnIndexWithinItsMatrix,lowerSubcolumnIndexWithinItsMatrix)];
					thisColumnCovariance.block(iUpperSubcolumn*L,iLowerSubcolumn*L,L,L) = subCovariance;
					thisColumnCovariance.block(iLowerSubcolumn*L,iUpperSubcolumn*L,L,L) = subCovariance.transpose(); // ...due to the symmetry of the covariance matrix
				}
			} // for(uint j=i;j<(d+1);j++)
		}
		
		columnsAutoCorrelationSum += thisColumnCovariance + predictedStackedChannelMatrix.col(iCol)*predictedStackedChannelMatrix.col(iCol).transpose();

	} // for(uint iCol=0;iCol<predictedStackedChannelMatrix.cols();iCol++)

	MatrixXd _Rx = noiseCovariance + _alphabetVariance*columnsAutoCorrelationSum;

    _filter = _Rx.inverse()*predictedStackedChannelMatrix*_alphabetVariance;

    VectorXd softEstimations = _filter.transpose()*observations;

    // required for nthSymbolVariance computing
    _channelMatrix = predictedStackedChannelMatrix;

    return softEstimations.segment(_detectionStart,_nSymbolsToBeDetected);
}

KalmanFilterAwareMMSEDetector* KalmanFilterAwareMMSEDetector::clone()
{
	return new KalmanFilterAwareMMSEDetector(*this);
}
