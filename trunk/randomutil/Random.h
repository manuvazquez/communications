#ifndef RANDOM_H
#define RANDOM_H 1

#include <stdint.h>
#include <complex>

using namespace std;

class Random
{

	protected:
		uint32_t _seed;

	public:
		Random (uint32_t seed) : _seed(seed) { };
		float randn();
		complex<float> complexRandn();
		int randab(int a, int b) { return (a+(int) ((double) b*rand_r(&_seed)/(RAND_MAX+1.0))); };
		double rand() { return ((double) rand_r(&_seed))/((double) RAND_MAX); }
		uint32_t getSeed() { return _seed; }
		void setSeed(uint32_t seed) { _seed = seed; }
};

#endif
