#ifndef RANDOM_H
#define RANDOM_H 1

#include <stdint.h>
#include <sys/time.h>
#include <complex>

class Random
{

	protected:
		uint32_t _seed;
		bool _havesmpl;
		double _smpl;

	public:
		Random () {struct timeval tv; gettimeofday(&tv, NULL); _seed = tv.tv_sec * tv.tv_usec;}
		Random (uint32_t seed) : _seed(seed),_havesmpl(false) {};
		double randn();

		std::complex<double> complexRandn();
		int randab(int a, int b) { return (a+(int) ((double) b*rand_r(&_seed)/(RAND_MAX+1.0))); };
		double rand() { return ((double) rand_r(&_seed))/((double) RAND_MAX); }
		uint32_t getSeed() { return _seed; }
		void setSeed(uint32_t seed) { _seed = seed; }
};

#endif
