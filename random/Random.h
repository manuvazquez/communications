#ifndef RANDOM_H
#define RANDOM_H

#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include <complex>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class Random
{

	protected:
		uint32_t _seed;
		bool _havesmpl;
		double _smpl;
		

// 		static void toOctaveFileStreamAux(const std::vector<std::vector <Random> >&matrix,std::string name,std::ofstream &f);
		static void toOctaveFileStreamAux3D(const std::vector<std::vector<std::vector <Random> > >&matrix,std::string name,std::ofstream &f);

	public:
		Random ():_havesmpl(false),_smpl(0.0)
		{
			struct timeval tv; 
			gettimeofday(&tv, NULL); 
			_seed = tv.tv_sec * tv.tv_usec;
		}
		Random (uint32_t seed) : _seed(seed),_havesmpl(false),_smpl(0.0) {};
		double randn();

		std::complex<double> complexRandn();
		int randab(int a, int b) { return (a+(int) ((double) (b-a+1)*rand_r(&_seed)/(RAND_MAX+1.0))); };
		double rand() { return ((double) rand_r(&_seed))/((double) RAND_MAX); }
		uint32_t getSeed() { return _seed; }
		void setSeed(uint32_t seed) { _seed = seed; }
		void setStoredSample(double sample) {_smpl = sample; _havesmpl = true; }
		
		static void toOctaveFileStream(const std::vector <Random> &vector,std::string name,std::ofstream &f);
		static void toOctaveFileStream(const std::vector<std::vector <Random> >&matrix,std::string name,std::ofstream &f);
		static void toOctaveFileStream(const std::vector<std::vector<std::vector <Random> > >&matrix,std::string name,std::ofstream &f);
};

#endif
