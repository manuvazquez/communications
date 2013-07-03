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


#ifndef KALMANFILTERAWAREMMSEDETECTOR_H
#define KALMANFILTERAWAREMMSEDETECTOR_H

#include <MMSEDetector.h>

#include<iostream>
#include <KalmanEstimator.h>

class KalmanFilterAwareMMSEDetector : public MMSEDetector
{
protected:
	KalmanEstimator const *_kalmanEstimator;
	const std::vector<double> _ARcoefficients;
	
	const bool _interferenceCancellation;
	
	class CovarianceId
	{
	protected:
		int _t1,_t2;
		uint _c1,_c2;
		
	public:
		
		CovarianceId(int t1,int t2,uint c1,uint c2):_t1(t1),_t2(t2),_c1(c1),_c2(c2)
		{
		}
		
		bool operator<(const CovarianceId &id) const
		{	
			if(_t1==id._t1)
				if(_t2==id._t2)
					if(_c1==id._c1)
						if(_c2==id._c2)
							return false; // they are equal
						else return (_c2<id._c2);
					else return (_c1<id._c1);
				else
					return (_t2<id._t2);
			else 
				return (_t1<id._t1);
		}
	};
	
public:
    KalmanFilterAwareMMSEDetector(uint rows, uint cols, double alphabetVariance,uint nSymbolsToBeDetected,KalmanEstimator *kalmanEstimator,std::vector<double> ARcoefficients,bool interferenceCancellation=false);
	
    virtual VectorXd detect(const VectorXd &observations, const MatrixXd &channelMatrix, const MatrixXd& noiseCovariance);
	virtual VectorXd detect2orderAndAboveARprocess(const VectorXd &observations, const MatrixXd &channelMatrix, const MatrixXd& noiseCovariance);
	
    virtual KalmanFilterAwareMMSEDetector* clone();
	
	void setKalmanEstimator(KalmanEstimator * kalmanEstimator) { _kalmanEstimator = kalmanEstimator; }
};

#endif // KALMANFILTERAWAREMMSEDETECTOR_H
