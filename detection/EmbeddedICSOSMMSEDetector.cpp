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


#include "EmbeddedICSOSMMSEDetector.h"

#include <math.h>
#include <TransmissionUtil.h>

// #define DEBUG

EmbeddedICSOSMMSEDetector::EmbeddedICSOSMMSEDetector(uint rows, uint cols, double alphabetVariance, uint nSymbolsToBeDetected, KalmanEstimator* kalmanEstimator, vector< double > ARcoefficients, bool sos)
:SOSMMSEDetector(rows,cols,alphabetVariance,nSymbolsToBeDetected,kalmanEstimator,ARcoefficients,false), _sos(sos)
{
}

VectorXd EmbeddedICSOSMMSEDetector::detect(const VectorXd &observations, const MatrixXd &channelMatrix, const MatrixXd& noiseCovariance)
{	
	assert(_ARcoefficients.size()==1);
	
	uint nRows = channelMatrix.rows();
	
	uint N = _kalmanEstimator->nInputs();
	int Nm = _kalmanEstimator->cols();
	uint m = uint(Nm)/N;
	uint L = _kalmanEstimator->rows();
	
	// index of the first stacked channel matrix column to be accounted for when computing the channel matrix covariance
	uint iFirstColumn  = N*(m-1);
	
	// smoothing lag "d" is inferred from the received channel matrix and the internal Kalman filter (whose number of rows is L(d+1))
	uint d = nRows/L -1;
	
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
	
	VectorXd previousSymbols = Util::toVector(_pastInterferingSymbols, columnwise);
	MatrixXd autoCorrelationSum = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
	
	for(uint iCol1=0;iCol1<iFirstColumn;iCol1++)
	{
		for(uint iCol2=0;iCol2<iFirstColumn;iCol2++)
		{
			
			MatrixXd thisColumnsCovariance = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
			
			for(uint iCol1Subcolumn=0;iCol1Subcolumn<=d;iCol1Subcolumn++)
			{
				int col1SubcolumnIndexWithinItsMatrix = iCol1 - (iCol1Subcolumn*N);
				
				if(col1SubcolumnIndexWithinItsMatrix<0)
					
					continue;
				
				for(uint iCol2Subcolumn=0;iCol2Subcolumn<=d;iCol2Subcolumn++)
				{
					int col2SubcolumnIndexWithinItsMatrix = iCol2 - (iCol2Subcolumn*N);
					
					if(col2SubcolumnIndexWithinItsMatrix<0)
					
						continue;
					
					MatrixXd subCovariance;
					
					if(iCol1Subcolumn<=iCol2Subcolumn)
					{
					
						subCovariance = Util::subMatrixFromVectorIndexes(predictedCovariances[iCol1Subcolumn],
												iCol2indexesWithinKFstateVector[col1SubcolumnIndexWithinItsMatrix],iCol2indexesWithinKFstateVector[col2SubcolumnIndexWithinItsMatrix])
												*pow(_ARcoefficients[0],double(iCol2Subcolumn-iCol1Subcolumn));						
					}else
					{
						subCovariance = Util::subMatrixFromVectorIndexes(predictedCovariances[iCol2Subcolumn],
												iCol2indexesWithinKFstateVector[col2SubcolumnIndexWithinItsMatrix],iCol2indexesWithinKFstateVector[col1SubcolumnIndexWithinItsMatrix]).transpose()
												*pow(_ARcoefficients[0],double(iCol1Subcolumn-iCol2Subcolumn));
					}
					
					thisColumnsCovariance.block(iCol1Subcolumn*L,iCol2Subcolumn*L,L,L) = subCovariance;
				}
			}
			
			if(_sos)
				autoCorrelationSum += (thisColumnsCovariance + predictedStackedChannelMatrix.col(iCol1)*predictedStackedChannelMatrix.col(iCol2).transpose())*previousSymbols(iCol1)*previousSymbols(iCol2);
			else
				autoCorrelationSum += (predictedStackedChannelMatrix.col(iCol1)*predictedStackedChannelMatrix.col(iCol2).transpose())*previousSymbols(iCol1)*previousSymbols(iCol2);
		}
	}
	
	for(uint iCol=iFirstColumn;iCol<predictedStackedChannelMatrix.cols();iCol++)
	{
		MatrixXd thisColumnCovariance = MatrixXd::Zero(predictedStackedChannelMatrix.rows(),predictedStackedChannelMatrix.rows());
		
		for(uint iUpperSubcolumn=0;iUpperSubcolumn<=d;iUpperSubcolumn++)
		{
			int upperSubcolumnIndexWithinItsMatrix = iCol - (iUpperSubcolumn*N);
			for(uint iLowerSubcolumn=iUpperSubcolumn;iLowerSubcolumn<=d;iLowerSubcolumn++)
			{
				int lowerSubcolumnIndexWithinItsMatrix = iCol - (iLowerSubcolumn*N);
				
				// correlation with a vector of zeros is zero
				if(upperSubcolumnIndexWithinItsMatrix>=Nm || upperSubcolumnIndexWithinItsMatrix<0 || lowerSubcolumnIndexWithinItsMatrix>=Nm  || lowerSubcolumnIndexWithinItsMatrix<0)
				{
					continue;
				}
				
				// only "predictedCovariances[0]" is needed because the difference between a column in a given time instant and the same column a number of time steps later is given by a factor (which is the power below)
				MatrixXd subCovariance = Util::subMatrixFromVectorIndexes(predictedCovariances[iUpperSubcolumn],
											iCol2indexesWithinKFstateVector[upperSubcolumnIndexWithinItsMatrix],iCol2indexesWithinKFstateVector[lowerSubcolumnIndexWithinItsMatrix])
											*pow(_ARcoefficients[0],double(iLowerSubcolumn-iUpperSubcolumn));
				
				thisColumnCovariance.block(iUpperSubcolumn*L,iLowerSubcolumn*L,L,L) = subCovariance;
				thisColumnCovariance.block(iLowerSubcolumn*L,iUpperSubcolumn*L,L,L) = subCovariance.transpose(); // ...due to the symmetry of the covariance matrix

			} // for(uint j=i;j<(d+1);j++)
		}
		
		
		if(_sos)
			columnsAutoCorrelationSum += thisColumnCovariance + predictedStackedChannelMatrix.col(iCol)*predictedStackedChannelMatrix.col(iCol).transpose();
		else
			columnsAutoCorrelationSum += predictedStackedChannelMatrix.col(iCol)*predictedStackedChannelMatrix.col(iCol).transpose();
		
	} // for(uint iCol=0;iCol<predictedStackedChannelMatrix.cols();iCol++)
	
	MatrixXd _Rx = noiseCovariance + _alphabetVariance*columnsAutoCorrelationSum + autoCorrelationSum;
	
    _filter = _Rx.inverse()*channelMatrix.block(0,N*(m-1),L*(d+1),channelMatrix.cols()-N*(m-1))*_alphabetVariance;
	
    VectorXd softEstimations = _filter.transpose()*observations;
	
    // required for nthSymbolVariance computing
    _channelMatrix = channelMatrix;

    return softEstimations.segment(_detectionStart,_nSymbolsToBeDetected);
}

EmbeddedICSOSMMSEDetector* EmbeddedICSOSMMSEDetector::clone()
{
	return new EmbeddedICSOSMMSEDetector(*this);
}
