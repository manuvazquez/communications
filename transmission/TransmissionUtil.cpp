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


#include "TransmissionUtil.h"

MatrixXd TransmissionUtil::channelMatrices2stackedChannelMatrix(std::vector< MatrixXd > matrices, uint m, uint start, uint d)
{
	// incorrect number of columns in the matrices
	assert((matrices[0].cols() % m)==0);

    uint nMatricesToStack = d - start + 1;

	// insufficient number of matrices
	assert(matrices.size()>= nMatricesToStack);
	
	uint nOutputs = matrices[0].rows();
	uint nInputs = matrices[0].cols()/m;

	MatrixXd res = MatrixXd::Zero(nOutputs*nMatricesToStack,nInputs*(m+nMatricesToStack-1));

    for(uint i=start,iStartingFromZero=0;i<=d;i++,iStartingFromZero++)
		res.block(iStartingFromZero*nOutputs,iStartingFromZero*nInputs,nOutputs,nInputs*m) = matrices[i];

    return res;
}