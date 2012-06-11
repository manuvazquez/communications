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

KalmanFilterAwareMMSEDetector::KalmanFilterAwareMMSEDetector(uint rows, uint cols, double alphabetVariance,uint nSymbolsToBeDetected,KalmanEstimator *kalmanEstimator,std::vector<double> ARcoefficients)
:MMSEDetector(rows,cols,alphabetVariance,nSymbolsToBeDetected),_kalmanEstimator(kalmanEstimator),_ARcoefficients(ARcoefficients)
{
	// this implementation is only meaningful for channels that evolve according with a first order AR process
	assert(ARcoefficients.size()==1);
}

VectorXd KalmanFilterAwareMMSEDetector::detect(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance)
{
// 	return detect2(observations,channelMatrix,noiseCovariance);
	
	uint nRows = channelMatrix.rows();
	
	uint N = _kalmanEstimator->nInputs();
	int Nm = _kalmanEstimator->cols();
	uint m = uint(Nm)/N;
	uint L = _kalmanEstimator->rows();
	
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
		// Regarding the noise variance, a value different from 0 (any one should do) must be passed or otherwise, the matrix inversion within the KF will give rise to NaN's
		MatrixXd foo = kalmanEstimatorClone->nextMatrix(VectorXd::Zero(L),MatrixXd::Zero(N,m),1.0);
		
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
	
	for(uint iCol=0;iCol<predictedStackedChannelMatrix.cols();iCol++)
	{
		MatrixXd thisColumnCovariance = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
		
		for(uint iUpperSubcolumn=0;iUpperSubcolumn<(d+1);iUpperSubcolumn++)
		{
			int upperSubcolumnOriginalColumn = iCol - (iUpperSubcolumn*N);
			for(uint iLowerSubcolumn=iUpperSubcolumn;iLowerSubcolumn<(d+1);iLowerSubcolumn++)
			{
				int lowerSubcolumnOriginalColumn = iCol - (iLowerSubcolumn*N);
				
				if(upperSubcolumnOriginalColumn>=Nm || upperSubcolumnOriginalColumn<0)
					continue;
				
				if(lowerSubcolumnOriginalColumn>=Nm || lowerSubcolumnOriginalColumn<0)
					continue;
				
				if(iUpperSubcolumn == iLowerSubcolumn) // => upperSubcolumnOriginalColumn == lowerSubcolumnOriginalColumn
					thisColumnCovariance.block(iUpperSubcolumn*L,iUpperSubcolumn*L,L,L) = Util::subMatrixFromVectorIndexes(predictedCovariances[iUpperSubcolumn],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn]);
				else
				{
// 					MatrixXd subCovariance = Util::subMatrixFromVectorIndexes(predictedCovariances[0],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn])
					MatrixXd subCovariance = Util::subMatrixFromVectorIndexes(predictedCovariances[0],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn],iCol2indexesWithinKFstateVector[lowerSubcolumnOriginalColumn])
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

VectorXd KalmanFilterAwareMMSEDetector::detect2(VectorXd observations, MatrixXd channelMatrix, const MatrixXd& noiseCovariance)
{	
	uint nRows = channelMatrix.rows();
	
	uint N = _kalmanEstimator->nInputs();
	uint Nm = _kalmanEstimator->cols();
	uint m = uint(Nm)/N;
	uint L = _kalmanEstimator->rows();
	
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
		// Regarding the noise variance, a value different from 0 (any one should do) must be passed or otherwise, the matrix inversion within the KF will give rise to NaN's
		kalmanEstimatorClone->nextMatrix(VectorXd::Zero(L),MatrixXd::Zero(N,m),1.0);
		
		predictedMatrices[i] = kalmanEstimatorClone->predictedMatrix();
		predictedCovariances[i] = kalmanEstimatorClone->getInternalPredictiveCovariance();
	}
	
	delete kalmanEstimatorClone;
	
	MatrixXd predictedStackedChannelMatrix = TransmissionUtil::channelMatrices2stackedChannelMatrix(predictedMatrices,m);

	MatrixXd columnsAutoCorrelationSum = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());

	// for every possible column index, the indexes within the KF state vector that give that column are obtained
	std::vector<std::vector<uint> > iCol2indexesWithinKFstateVector(Nm);
	
	for(uint iCol=0;iCol<Nm;iCol++)
		iCol2indexesWithinKFstateVector[iCol] = _kalmanEstimator->colIndexToIndexesWithinKFstateVector(iCol,_ARcoefficients.size()-1);
	
	// ------------------------
	
	std::map<CovarianceId,MatrixXd> covariancesMap;
	
#ifdef DEBUG
	cout << "predictedCovariances[0]: " << endl << predictedCovariances[0] << endl;
#endif


	for(int t1=0;uint(t1)<_ARcoefficients.size();t1++)
		for(int t2=0;uint(t2)<_ARcoefficients.size();t2++)
			for(uint i=0;i<Nm;i++)
				for(uint j=0;j<i;j++)
// 				for(uint j=i+1;j<Nm;j++)
				{
					covariancesMap[CovarianceId(-t1,-t2,i,j)]= Util::subMatrixFromVectorIndexes(predictedCovariances[0],
													_kalmanEstimator->colIndexToIndexesWithinKFstateVector(i,_ARcoefficients.size()-1-t1),
													_kalmanEstimator->colIndexToIndexesWithinKFstateVector(j,_ARcoefficients.size()-1-t2));
#ifdef DEBUG
					cout << "adding covariance for t1 = " << -t1 << ",t2 = " << -t2 << ",i = " << i << ",j = " << j << endl;
					cout << covariancesMap[CovarianceId(-t1,-t2,i,j)] << endl;
#endif
				}
				
#ifdef DEBUG
	cout << "second step ========================= " << endl;
#endif
	
	for(int t1 = -(_ARcoefficients.size()-1);t1<=0;t1++)
		for(int t2=1;t2<=int(d);t2++)
			for(uint i=0;i<Nm;i++)
				for(uint j=0;j<i;j++)
// 				for(uint j=i+1;j<Nm;j++)
				{
					MatrixXd covariance = MatrixXd::Zero(L,L);
					for(int r=1;r<=int(_ARcoefficients.size());r++)
					{
						assert(covariancesMap.find(CovarianceId(t1,t2-r,i,j))!=covariancesMap.end());
						covariance += _ARcoefficients[r-1]*covariancesMap[CovarianceId(t1,t2-r,i,j)];
#ifdef DEBUG
						cout << "adding covariance for t1 = " << t1 << ",t2 = " << t2 << ",i = " << i << ",j = " << j << endl;
						cout << "multiplying " << _ARcoefficients[r-1] << " times " << endl << covariancesMap[CovarianceId(t1,t2-r,i,j)] << endl << "gives" << endl << covariance << endl;
#endif
					}
					covariancesMap[CovarianceId(t1,t2,i,j)] = covariance;
				}

#ifdef DEBUG
	cout << "3rd step ========================= " << endl;
#endif
	
	for(int t1=1;t1<=int(d);t1++)
		for(int t2=t1+1;t2<=int(d);t2++)
			for(uint i=0;i<Nm;i++)
				for(uint j=0;j<i;j++)
// 				for(uint j=i+1;j<Nm;j++)
				{
					MatrixXd covariance = MatrixXd::Zero(L,L);
					for(uint r=1;r<=_ARcoefficients.size();r++)
					{
						assert(covariancesMap.find(CovarianceId(t1-int(r),t2,i,j))!=covariancesMap.end());
						covariance += _ARcoefficients[r-1]*covariancesMap[CovarianceId(t1-r,t2,i,j)];
#ifdef DEBUG
						cout << "adding covariance for t1 = " << t1 << ",t2 = " << t2 << ",i = " << i << ",j = " << j << endl;
						cout << "multiplying " << _ARcoefficients[r-1] << " times " << endl << covariancesMap[CovarianceId(t1-r,t2,i,j)] << endl << "gives" << endl << covariance << endl;
#endif
					}
					covariancesMap[CovarianceId(t1,t2,i,j)] = covariance;
				}

#ifdef DEBUG
// 	cout << "printing everything================" << endl;
// 	for ( std::map<CovarianceId,MatrixXd>::const_iterator iter = covariancesMap.begin();
// 	iter != covariancesMap.end(); ++iter )
// 	{
// 		cout << iter->first._t1 << " " << iter->first._t2 << " " << iter->first._c1 << " " << iter->first._c2 << endl << iter->second << '\n';
// 	}
// 	cout << endl << "================" << endl;
	
	cout << "size of the map is: " << covariancesMap.size() << endl;
	
	getchar();
#endif
	
	// -------------
	
	for(uint iCol=0;iCol<predictedStackedChannelMatrix.cols();iCol++)
	{
		MatrixXd thisColumnCovariance = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
		
		for(uint iUpperSubcolumn=0;iUpperSubcolumn<(d+1);iUpperSubcolumn++)
		{
			int upperSubcolumnOriginalColumn = iCol - (iUpperSubcolumn*N);
			for(uint iLowerSubcolumn=iUpperSubcolumn;iLowerSubcolumn<(d+1);iLowerSubcolumn++)
			{
				int lowerSubcolumnOriginalColumn = iCol - (iLowerSubcolumn*N);
				
				if(upperSubcolumnOriginalColumn>=int(Nm) || upperSubcolumnOriginalColumn<0)
					continue;
				
				if(lowerSubcolumnOriginalColumn>=int(Nm) || lowerSubcolumnOriginalColumn<0)
					continue;
				
				if(iUpperSubcolumn == iLowerSubcolumn) // => upperSubcolumnOriginalColumn == lowerSubcolumnOriginalColumn
					thisColumnCovariance.block(iUpperSubcolumn*L,iUpperSubcolumn*L,L,L) = Util::subMatrixFromVectorIndexes(predictedCovariances[iUpperSubcolumn],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn]);
				else
				{
// 					MatrixXd subCovariance = Util::subMatrixFromVectorIndexes(predictedCovariances[0],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn])
// 					MatrixXd subCovariance = Util::subMatrixFromVectorIndexes(predictedCovariances[0],iCol2indexesWithinKFstateVector[upperSubcolumnOriginalColumn],iCol2indexesWithinKFstateVector[lowerSubcolumnOriginalColumn])
// 												*pow(_ARcoefficients[0],double(iUpperSubcolumn+iLowerSubcolumn));
					MatrixXd subCovariance = covariancesMap[CovarianceId(iUpperSubcolumn,iLowerSubcolumn,upperSubcolumnOriginalColumn,lowerSubcolumnOriginalColumn)];
					thisColumnCovariance.block(iUpperSubcolumn*L,iLowerSubcolumn*L,L,L) = subCovariance;
					thisColumnCovariance.block(iLowerSubcolumn*L,iUpperSubcolumn*L,L,L) = subCovariance.transpose(); // ...due to the symmetry of the covariance matrix
#ifdef DEBUG
					cout << "iUpperSubcolumn = " << iUpperSubcolumn << " iLowerSubcolumn = " << iLowerSubcolumn << " upperSubcolumnOriginalColumn = " << upperSubcolumnOriginalColumn << " lowerSubcolumnOriginalColumn " << lowerSubcolumnOriginalColumn << endl;
					cout << "traditional" << endl << subCovariance << endl;
					cout << "new" << endl << covariancesMap[CovarianceId(iUpperSubcolumn,iLowerSubcolumn,upperSubcolumnOriginalColumn,lowerSubcolumnOriginalColumn)] << endl;
#endif
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

KalmanFilterAwareMMSEDetector* KalmanFilterAwareMMSEDetector::clone()
{
	return new KalmanFilterAwareMMSEDetector(*this);
}

// MatrixXd KalmanFilterAwareMMSEDetector::subCovar(const MatrixXd &covariance,int t1,int t2,uint i,uint j) const
// {
// 	assert(t1>-int(_ARcoefficients.size()));
// 	assert(t2>-int(_ARcoefficients.size()));
// 	
// 	if(t1<=0 && t2<=0)
// 		return Util::subMatrixFromVectorIndexes(covariance,
// 										 _kalmanEstimator->colIndexToIndexesWithinKFstateVector(i,_ARcoefficients.size()-1+t1),
// 										 _kalmanEstimator->colIndexToIndexesWithinKFstateVector(j,_ARcoefficients.size()-1+t2));
// 	else
// 	{
// 	}
// }