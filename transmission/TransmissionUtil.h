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


#ifndef TRANSMISSIONUTIL_H
#define TRANSMISSIONUTIL_H

#include <types.h>
#include <vector>

class TransmissionUtil
{
public:
    static MatrixXd channelMatrices2stackedChannelMatrix(std::vector<MatrixXd> matrices,uint m,uint start,uint d);
    static MatrixXd channelMatrices2stackedChannelMatrix(std::vector<MatrixXd> matrices,uint m)
    {
        return channelMatrices2stackedChannelMatrix(matrices,m,0,matrices.size()-1);
    }
    static MatrixXd channelMatrices2stackedChannelMatrix(std::vector<MatrixXd> matrices,uint m,uint d)
    {
        return channelMatrices2stackedChannelMatrix(matrices,m,0,d);
    }
};

#endif // TRANSMISSIONUTIL_H
