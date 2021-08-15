#ifndef RANDOMBITS_H_
#define RANDOMBITS_H_

#include <string>
#include <random>
#include <vector>
#include <gmpxx.h>

using std::string;
using std::mt19937;
using std::mt19937_64;
using std::vector;

class RandomBits {
	private:
		string binary(unsigned x, uint32_t length);
		mt19937 randomEngine{};
		mt19937_64 randomEngine2{};
		gmp_randclass g;
		void seedEngine(std::random_device&);
		void seedEngine2(std::random_device&);
	public:
		void seedEngines();
		RandomBits();
		string GenerateRandomBits(uint32_t size);
		bool generateWeightedRandomBit(long double posWt, long double totWt);
		//bool generateMPFWeightedRandomBit(mpf_class& posWt, mpf_class& totWt);
		bool generateWeightedRandomBit(mpf_class& posWt, mpf_class& totWt);
		uint64_t getRandInt(std::uniform_int_distribution<uint64_t>& uid);
		double_t getRandReal(std::uniform_real_distribution<double_t>& urd);
		long double getRandReal(std::uniform_real_distribution<long double>& urd);
};

#endif
