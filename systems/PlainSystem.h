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
#include <LinearFilterAwareNoiseVarianceAdjustingKalmanEstimator.h>

#include <KalmanFilterAwareMMSEDetector.h>

#include <KnownChannelChannelMatrixEstimator.h>

#include <MMSEDetector.h>

#include <FlatPowerProfile.h>

#include <LinearFilterBasedAlgorithm.h>
#include <LinearFilterBasedAlgorithmWithKnownChannel.h>
#include <LinkedKalmanFilterAndLinearFilterBasedAlgorithm.h>
#include <KnownSymbolsKalmanBasedChannelEstimatorAlgorithm.h>
#include <KalmanFilterAwareMMSEBasedAlgorithm.h>


class PlainSystem : public BaseSystem
{

protected:
	KalmanEstimator *_kalmanEstimator;
	KnownSymbolsKalmanEstimator *_knownSymbolsKalmanEstimator;
	
	MMSEDetector *_MMSEdetector;
	KnownChannelChannelMatrixEstimator *_knownChannelChannelMatrixEstimator;
	
	KalmanFilterAwareMMSEDetector *_kalmanFilterAwareMMSEDetector;
	
	
    virtual void addAlgorithms();

public:
    PlainSystem();
    virtual ~PlainSystem();
};

#endif // PLAINSYSTEM_H
