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


#ifndef PLAINSYSTEM_H
#define PLAINSYSTEM_H

#include <BaseSystem.h>

#include <KalmanEstimator.h>
#include <KnownSymbolsKalmanEstimator.h>

#include <SOSMMSEDetector.h>

#include <KnownChannelChannelMatrixEstimator.h>

#include <MMSEDetector.h>

#include <FlatPowerProfile.h>

#include <LinearFilterKFBasedAlgorithm.h>
#include <KnownSymbolsKalmanBasedChannelEstimatorAlgorithm.h>
#include <SOSMMSEBasedAlgorithm.h>


class PlainSystem : public BaseSystem
{

protected:
	KalmanEstimator *_kalmanEstimator;
	KnownSymbolsKalmanEstimator *_knownSymbolsKalmanEstimator;
	
	/**
	 * @brief interference-cancelating MMSE detector
	 **/
	MMSEDetector *_ICMMSEdetector;
	
	KnownChannelChannelMatrixEstimator *_knownChannelChannelMatrixEstimator;
	
	/**
	 * @brief SOS-MMSE detector
	 **/
	SOSMMSEDetector *_SOSMMSEDetector;
	
	/**
	 * @brief interference-cancelating SOS-MMSE detector (an MMSE detector taking advantage of second-order statistics meant to be used after interference cancellation)
	 **/
	SOSMMSEDetector *_ICSOSMMSEDetector;
	
    virtual void addAlgorithms();

public:
    PlainSystem();
    virtual ~PlainSystem();
};

#endif // PLAINSYSTEM_H
