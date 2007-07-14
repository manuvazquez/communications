#ifndef PARAMETERS_DEFINED
vector<int> __MMSEforwardSmoothing;
__MMSEforwardSmoothing.push_back(d);__MMSEforwardSmoothing.push_back(5);__MMSEforwardSmoothing.push_back(10);
__MMSEforwardSmoothing.push_back(20);__MMSEforwardSmoothing.push_back(50);

vector<LinearDetector *> __MMSEforwardSmoothingDetectors;

for(uint iSmoothing=0;iSmoothing<__MMSEforwardSmoothing.size();iSmoothing++)
	__MMSEforwardSmoothingDetectors.push_back(new MMSEDetector(L*(c+__MMSEforwardSmoothing[iSmoothing]+1),N*(__MMSEforwardSmoothing[iSmoothing]+1),pam2.Variance(),N*(d+1),0));
#else
for(uint iSmoothing=0;iSmoothing<__MMSEforwardSmoothing.size();iSmoothing++)
{
	char buffer[SPRINTF_BUFFER];

	sprintf(buffer," e = %d",__MMSEforwardSmoothing[iSmoothing]);

	algorithms.push_back(new LinearFilterBasedMKFAlgorithm("MKF (MMSE)" + string(buffer),pam2,L,N,lastSymbolVectorInstant,m,&kalmanEstimator,__MMSEforwardSmoothingDetectors[iSmoothing],preamble,c,d,__MMSEforwardSmoothing[iSmoothing],nParticles,&algoritmoRemuestreo,powerProfile.Means(),powerProfile.Variances(),ARcoefficients[0],firstSampledChannelMatrixVariance,ARvariance));
}
#endif
